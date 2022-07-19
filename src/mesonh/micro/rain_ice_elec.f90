!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #########################
       MODULE MODI_RAIN_ICE_ELEC
!      #########################
!
INTERFACE
      SUBROUTINE RAIN_ICE_ELEC (OSEDIC, HSUBG_AUCV, OWARM,                            &
                                KSPLITR, PTSTEP, KMI, KRR,                            &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                                PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                                PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                                PINPRS, PINPRG, PSIGS,                                &
                                PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,           &
                                PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,           &
                                PSEA, PTOWN,                                          &
                                PRHT, PRHS, PINPRH, PQHT, PQHS                        )
!
!
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                                   ! integration for rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRC    ! Cloud instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRR    ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D  ! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D   ! Rain evap profile
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRS    ! Snow instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRG    ! Graupel instant precip
!
! Charge Mixing Ratio (CMR) (C/kg)
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQPIT   ! Positive ion (Nb/kg) at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQNIT   ! Negative ion (Nb/kg) at t
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQCT    ! Cloud water CMR at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQRT    ! Rain water CMR at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQIT    ! Pristine ice CMR at t
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQST    ! Snow/aggregate CMR at t
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQGT    ! Graupel CMR at t
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQPIS   ! Positive ion source 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQNIS   ! Negative ion source 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQCS    ! Cloud water CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQRS    ! Rain water CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQIS    ! Pristine ice CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQSS    ! Snow/aggregate CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQGS    ! Graupel CMR source
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PSEA
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PTOWN
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(INOUT) :: PINPRH  ! Hail instant precip
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PQHT    ! Hail CMR at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PQHS    ! Hail CMR source
!
END SUBROUTINE RAIN_ICE_ELEC
END INTERFACE
END MODULE MODI_RAIN_ICE_ELEC
!
!     ######spl
      SUBROUTINE RAIN_ICE_ELEC (OSEDIC, HSUBG_AUCV, OWARM,                            &
                                KSPLITR, PTSTEP, KMI, KRR,                            &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                                PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                                PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                                PINPRS, PINPRG, PSIGS,                                &
                                PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,           &
                                PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,           &
                                PSEA, PTOWN,                                          &
                                PRHT, PRHS, PINPRH, PQHT, PQHS                        )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources
!!           and the cloud electrification
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the slow microphysical sources
!!    which can be computed explicitly
!!
!!
!!**  METHOD
!!    ------
!!      The autoconversion computation follows Kessler (1969).
!!      The sedimentation rate is computed with a time spliting technique and 
!!    an upstream scheme, written as a difference of non-advective fluxes. This
!!    source term is added to the future instant ( split-implicit process ).
!!      The others microphysical processes are evaluated at the central instant 
!!    (split-explicit process ): autoconversion, accretion and rain evaporation.
!!      These last 3 terms are bounded in order not to create negative values 
!!    for the water species at the future instant.
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
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XCI                ! Ci (solid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                               function over liquid water
!!          XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure
!!                               function over solid ice
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!
!!    AUTHOR
!!    ------
!!      C. Barthe, G. Molinie, J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2002
!!      Modifications
!!         C. Barthe (LACy)  Nov. 2009 : update to V4.8.1
!!         M. Chong      26/01/10  Add Small ions parameters
!!         J-P Pinty     31/03/11  Add hail
!!         C. Lac        2011 : Adaptation to FIT temporal scheme
!!         B. Tsenova    June 2012 Add new NI parameterizations
!!         C. Barthe     June 2012 Dependance of RAR on the RELATIVE terminal velocity
!!         M. Chong      06/08/13  Add "Beard" effect (ELEC=>MICROPHYSICS) 
!!         J-P Pinty     21/08/13  Correction of the process limitation algo.
!!                                 SIGN(MIN(ABS ...
!!                                 Correction in elec_update_qd
!!                                 Correction of hail charge transfer
!!                                 Add hail growth charging processes
!!         J-P Pinty     26/08/13  Add "Beard" effect control (ELEC=>MICROPHYS)
!!                                 for sedimentation
!!         J-P Pinty     26/09/13  Add tabulated treatment of SAUN1 and SAUN2
!!         J-P Pinty     30/09/13  Remove call to MOMG function
!!         J-P Pinty     25/10/13  Add "Latham" effect for aggregation process
!!         M. Chong      31/10/13  Add other tabulated treatment and recode
!!         M. Chong      15/11/13  Bug in the computation of RGWETH (wrong sign)
!!         J-P Pinty     25/04/14  Many bugs with ZWQ1(:,...) = 0.0
!!         J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!         J.Escobar : 10/2017 : for real*4 , limit exp() in RAIN_ICE_ELEC_SLOW with XMNH_HUGE_12_LOG
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  P .Wautelet 09/03/2020: add missing budgets for electricity
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbu_enable,                                                 &
                                lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri, &
                                lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,             &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, &
                                NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1,            &
                                tbudgets
USE MODD_CONF
USE MODD_CST
USE MODD_ELEC_DESCR
USE MODD_ELEC_n
USE MODD_ELEC_PARAM
USE MODD_LES
USE MODE_ll
USE MODD_NSV,             ONLY: NSV_ELECBEG, NSV_ELECEND ! Scalar variables for budgets
USE MODD_PARAMETERS
USE MODD_PARAM_ICE
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_REF,             ONLY: XTHVREFZ

use mode_budget,          only: Budget_store_add, Budget_store_init, Budget_store_end
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
use mode_tools,           only: Countjv

USE MODI_MOMG

IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),       INTENT(IN)    :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                 !   form by warm processes
                                                 !      (Kessler scheme)
!
INTEGER,                INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                   INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                INTENT(IN)    :: KMI     ! Model index 
INTEGER,                INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCLDFR! Convective Mass Flux Cloud fraction
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRC   ! Cloud instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D ! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! Rain evap profile
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRG   ! Graupel instant precip
!
! Charge Mixing Ratio (CMR) (C/kg)
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQPIT   ! Positive ion (Nb/kg) at t  
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQNIT   ! Negative ion (Nb/kg) at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQCT    ! Cloud water CMR at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQRT    ! Rain water CMR at t 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQIT    ! Pristine ice CMR at t
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQST    ! Snow/aggregate CMR at t
REAL, DIMENSION(:,:,:), INTENT(IN) :: PQGT    ! Graupel CMR at t
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQPIS   ! Positive ion source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQNIS   ! Negative ion source 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQCS    ! Cloud water CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQRS    ! Rain water CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQIS    ! Pristine ice CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQSS    ! Snow/aggregate CMR source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQGS    ! Graupel CMR source
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PSEA
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PTOWN
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(INOUT) :: PINPRH  ! Hail instant precip
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PQHT    ! Hail CMR at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PQHS    ! Hail CMR source
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: JJ            ! Loop index for the interpolation
INTEGER :: JI            ! Loop index for the interpolation
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
!
INTEGER :: ISEDIMR,ISEDIMC, ISEDIMI, ISEDIMS, ISEDIMG, ISEDIMH, &
           INEGT, IMICRO ! Case number of sedimentation, T>0 (for HEN)
                         ! and r_x>0 locations
INTEGER :: IGRIM, IGACC, IGDRY ! Case number of riming, accretion and dry growth
                               ! locations
INTEGER :: IGWET, IHAIL   ! wet growth locations and case number
!
LOGICAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
    :: GSEDIMR, GSEDIMC, GSEDIMI, GSEDIMS, GSEDIMG, GSEDIMH ! Test where to compute the SED processes
LOGICAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                          :: GNEGT  ! Test where to compute the HEN process
LOGICAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                          :: GMICRO ! Test where to compute all processes
LOGICAL, DIMENSION(:), ALLOCATABLE :: GRIM ! Test where to compute riming
LOGICAL, DIMENSION(:), ALLOCATABLE :: GACC ! Test where to compute accretion
LOGICAL, DIMENSION(:), ALLOCATABLE :: GDRY ! Test where to compute dry growth
LOGICAL, DIMENSION(:), ALLOCATABLE :: GWET  ! Test where to compute wet growth
LOGICAL, DIMENSION(:), ALLOCATABLE :: GHAIL ! Test where to compute hail growth
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2       ! Vectors of indices for
                                                        ! interpolations
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for 
                                                        ! interpolations
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZW  ! work array
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
      :: ZPRCS, ZPRRS, ZPRSS, ZPRGS, ZPRHS   ! Mixing ratios created during the time step
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZWSED         ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZWSEDW1       ! sedimentation speed
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZWSEDW2       ! sedimentation speed
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2))                   &
                                  :: ZCONC_TMP     ! Weighted concentration
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZT            ! Temperature
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::  &
                                     ZRAY,   & ! Cloud Mean radius
                                     ZLBC,   & ! XLBC weighted by sea fraction
                                     ZFSEDC
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
      :: ZPQRS, ZPQSS, ZPQGS, ZPQHS   ! Charge Mixing ratios created during the time step
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZWSEDQ         ! sedimentation fluxes for charge
REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRGS    ! Graupel m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRHS    ! Hail m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZTHS    ! Theta source
REAL, DIMENSION(:), ALLOCATABLE :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
!
REAL, DIMENSION(:), ALLOCATABLE &
               :: ZRHODREF, & ! RHO Dry REFerence
                  ZRHODREFC,& ! RHO Dry REFerence
                  ZRHODREFR,& ! RHO Dry REFerence
                  ZRHODREFI,& ! RHO Dry REFerence
                  ZRHODREFS,& ! RHO Dry REFerence
                  ZRHODREFG,& ! RHO Dry REFerence
                  ZRHODREFH,& ! RHO Dry REFerence
                  ZRHODJ,   & ! RHO times Jacobian
                  ZZT,      & ! Temperature
                  ZPRES,    & ! Pressure
                  ZEXNREF,  & ! EXNer Pressure REFerence
                  ZZW,      & ! Work array
                  ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                  ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                  ZUSW,     & ! Undersaturation over water
                  ZSSI,     & ! Supersaturation over ice
                  ZLBDAI,   & ! Slope parameter of the pristine ice distribution
                  ZLBDAR,   & ! Slope parameter of the raindrop  distribution
                  ZLBDAS,   & ! Slope parameter of the aggregate distribution
                  ZLBDAG,   & ! Slope parameter of the graupel   distribution
                  ZLBDAH,   & ! Slope parameter of the hail      distribution
                  ZRDRYG,   & ! Dry growth rate of the graupeln
                  ZRWETG,   & ! Wet growth rate of the graupeln
                  ZAI,      & ! Thermodynamical function
                  ZCJ,      & ! Function to compute the ventilation coefficient
                  ZKA,      & ! Thermal conductivity of the air
                  ZDV,      & ! Diffusivity of water vapor in the air
                  ZSIGMA_RC,& ! Standard deviation of rc at time t  
                  ZCF,      & ! Cloud fraction
                  ZCC,      & ! terminal velocity
                  ZFSEDC1D, & ! For cloud sedimentation
                  ZWLBDC,   & ! Slope parameter of the droplet distribution
                  ZCONC,    & ! Concentration des aérosols
                  ZRAY1D,   & ! Mean radius
                  ZWLBDA      ! Libre parcours moyen
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW1 ! Work arrays
REAL                              :: ZTIMAUTIC
REAL, DIMENSION(SIZE(XRTMIN))     :: ZRTMIN
!
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
LOGICAL, DIMENSION(:,:),ALLOCATABLE :: GELEC ! Logical of work for elec
REAL, DIMENSION(:),   ALLOCATABLE :: ZRSMIN_ELEC  ! Limit value of ZRXS where charge is available
REAL, DIMENSION(:),   ALLOCATABLE :: ZVECQ4, &  ! Work
                                     ZVECQ5, &  ! vectors for
                                     ZVECQ6, &  ! interpolations
                                     ZVECQ7     ! (electrification)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWQ1      ! Work array for electrification
REAL, DIMENSION(:),   ALLOCATABLE :: ZWQ3,ZWQ4 ! Work arrays for electrification

REAL, DIMENSION(:), ALLOCATABLE :: ZQPIT  ! Positive ion (kg^-1) at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZQNIT  ! Negative ion (kg^-1) at t
REAL, DIMENSION(:), ALLOCATABLE :: ZQCT   ! Cloud water CMR at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZQRT   ! Rain water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZQIT   ! Pristine ice CMR at t
REAL, DIMENSION(:), ALLOCATABLE :: ZQST   ! Snow/aggregate CMR at t
REAL, DIMENSION(:), ALLOCATABLE :: ZQGT   ! Graupel CMR at t
REAL, DIMENSION(:), ALLOCATABLE :: ZQHT   ! Hail CMR at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZQPIS  ! Positive ion source 
REAL, DIMENSION(:), ALLOCATABLE :: ZQNIS  ! Negative ion source 
REAL, DIMENSION(:), ALLOCATABLE :: ZQCS   ! Cloud water CMR source
REAL, DIMENSION(:), ALLOCATABLE :: ZQRS   ! Rain water CMR source
REAL, DIMENSION(:), ALLOCATABLE :: ZQIS   ! Pristine ice CMR source
REAL, DIMENSION(:), ALLOCATABLE :: ZQSS   ! Snow/aggregate CMR source
REAL, DIMENSION(:), ALLOCATABLE :: ZQGS   ! Graupel CMR source
REAL, DIMENSION(:), ALLOCATABLE :: ZQHS   ! Hail CMR source
!
! Charge diameter relation 
REAL, DIMENSION(:), ALLOCATABLE :: ZECT   ! Cloud water at t
REAL, DIMENSION(:), ALLOCATABLE :: ZERT   ! Rain water at t
REAL, DIMENSION(:), ALLOCATABLE :: ZEIT   ! Pristine ice at t
REAL, DIMENSION(:), ALLOCATABLE :: ZEST   ! Snow/aggregate at t
REAL, DIMENSION(:), ALLOCATABLE :: ZEGT   ! Graupel at t
REAL, DIMENSION(:), ALLOCATABLE :: ZEHT   ! Hail at t
! 
REAL, DIMENSION(:), ALLOCATABLE :: ZECS   ! Cloud water at t+dt 
REAL, DIMENSION(:), ALLOCATABLE :: ZERS   ! Rain water at t+dt
REAL, DIMENSION(:), ALLOCATABLE :: ZEIS   ! Pristine ice at t+dt
REAL, DIMENSION(:), ALLOCATABLE :: ZESS   ! Snow/aggregate at t+dt
REAL, DIMENSION(:), ALLOCATABLE :: ZEGS   ! Graupel at t+dt
REAL, DIMENSION(:), ALLOCATABLE :: ZEHS   ! Hail at t+dt
!
REAL, DIMENSION(:), ALLOCATABLE :: ZDELTALWC ! Gap between LWC and a critical LWC
REAL, DIMENSION(:), ALLOCATABLE :: ZLWCC     ! Critical LWC in NI charging
REAL, DIMENSION(:), ALLOCATABLE :: ZFT       ! Fct depending on temperature
!
! Non-inductive charging process following Saunders et al. (1991) / EW
REAL, DIMENSION(:), ALLOCATABLE :: ZEW      ! Effective liquid water content
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK  ! constant B _______________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM  !  d_i exponent ____________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN  !  v_g/s-v_i________________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM  !  d_s exponent ____________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN  !  v_g-v_s _________________________
REAL, DIMENSION(:), ALLOCATABLE :: ZFQIAGGS, ZFQIDRYGBS
REAL, DIMENSION(:), ALLOCATABLE :: ZLBQSDRYGB1S, ZLBQSDRYGB2S, ZLBQSDRYGB3S
!
! Non-inductive charging process following Saunders and Peck (1998) / RAR
REAL, DIMENSION(:), ALLOCATABLE :: ZVGMEAN  ! Mean velocity of graupel
REAL, DIMENSION(:), ALLOCATABLE :: ZVSMEAN  ! Mean velocity of snow
REAL, DIMENSION(:), ALLOCATABLE :: ZRHOCOR  ! Density correction for fallspeed
REAL, DIMENSION(:), ALLOCATABLE :: ZRAR     ! Rime accretion rate
REAL, DIMENSION(:), ALLOCATABLE :: ZRAR_CRIT  ! Critical RAR
REAL, DIMENSION(:), ALLOCATABLE :: ZDQRAR_IS   ! q= f(RAR,T) in Saunders and Peck's equation
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IS  !  d_i exponent ____________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IS  !  v_g/s-v_i________________________
REAL, DIMENSION(:), ALLOCATABLE :: ZDQRAR_IG   ! q= f(RAR,T) in Saunders and Peck's equation
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IG  !  d_i exponent ____________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IG  !  v_g/s-v_i________________________
REAL, DIMENSION(:), ALLOCATABLE :: ZDQRAR_SG   ! q= f(RAR,T) in Saunders and Peck's equation
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK_SG  ! constant B _______________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM_SG  !  d_s exponent ____________________
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN_SG  !  v_g-v_s _________________________
!
! Non-inductive charging process following Takahashi (1978)
INTEGER :: IGTAKA ! Case number of charge separation for Takahashi param.
LOGICAL, DIMENSION(:), ALLOCATABLE :: GTAKA      ! Test where to compute charge
                                                 ! separation for Takahashi param.
REAL,    DIMENSION(:), ALLOCATABLE :: ZDQTAKA_OPT ! Optimized array of separated charge
!
INTEGER :: IGSAUN ! Case number of charge separation for Saunders param.
LOGICAL, DIMENSION(:), ALLOCATABLE :: GSAUN      ! Test where to compute charge
                                                 ! separation for Saunders param.
REAL,    DIMENSION(:), ALLOCATABLE :: ZDQLWC_OPT ! Optimized array of separated charge
REAL,    DIMENSION(:), ALLOCATABLE :: ZDQLWC     ! q=f(LWC,T)
!
! Inductive charging process (Ziegler et al., 1991)
INTEGER                            :: IIND      ! Case number of inductive process
LOGICAL, DIMENSION(:), ALLOCATABLE :: GIND      ! Test where to compute inductive process
REAL,    DIMENSION(:), ALLOCATABLE :: ZRATE_IND ! Charge transfer rate during inductive process
REAL,    DIMENSION(:), ALLOCATABLE :: ZEFIELDW  ! Vertical component of the electric field
!
! Latham's effect
REAL,    DIMENSION(:), ALLOCATABLE :: ZLATHAMIAGGS ! E Function to simulate
                                                   ! enhancement of IAGGS
REAL,    DIMENSION(:), ALLOCATABLE :: ZEFIELDU  ! Horiz. component of the electric field
REAL,    DIMENSION(:), ALLOCATABLE :: ZEFIELDV  ! Horiz. component of the electric field
!
REAL, DIMENSION(:), ALLOCATABLE :: ZLIMIT, ZAUX, ZAUX1
REAL, DIMENSION(:), ALLOCATABLE :: ZCOLIS   ! Collection efficiency between ice and snow
REAL, DIMENSION(:), ALLOCATABLE :: ZCOLIG   ! Collection efficiency between ice and graupeln
REAL, DIMENSION(:), ALLOCATABLE :: ZCOLSG   ! Collection efficiency between snow and graupeln
REAL :: ZRHO00, ZCOR00    ! Surface reference air density
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PZZ,3) - JPVEXT
!
ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
ZCOR00 = ZRHO00**XCEXVT
!
!
!*       2.     COMPUTES THE SLOW COLD PROCESS SOURCES
!   	        --------------------------------------
!
!*       2.1    compute the ice nucleation
!
CALL RAIN_ICE_ELEC_NUCLEATION
!
!
!*       2.2    allocations
!
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
GMICRO(:,:,:) = .FALSE.

IF ( KRR == 7 ) THEN
  GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                PRCT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(2) .OR. &
                PRRT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(3) .OR. &
                PRIT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(4) .OR. &
                PRST(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(5) .OR. &
                PRGT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(6) .OR. &
                PRHT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(7)
ELSE IF( KRR == 6 ) THEN
  GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                PRCT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(2) .OR. &
                PRRT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(3) .OR. &
                PRIT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(4) .OR. &
                PRST(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(5) .OR. &
                PRGT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(6)
END IF

IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
IF (IMICRO > 0) THEN
  ALLOCATE(ZRVT(IMICRO)) 
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))
  ALLOCATE(ZRIT(IMICRO))
  ALLOCATE(ZRST(IMICRO))
  ALLOCATE(ZRGT(IMICRO))
 IF (KRR == 7) ALLOCATE(ZRHT(IMICRO))
  ALLOCATE(ZCIT(IMICRO))
  ALLOCATE(ZRVS(IMICRO))
  ALLOCATE(ZRCS(IMICRO))
  ALLOCATE(ZRRS(IMICRO))
  ALLOCATE(ZRIS(IMICRO))
  ALLOCATE(ZRSS(IMICRO))
  ALLOCATE(ZRGS(IMICRO))
 IF (KRR == 7) ALLOCATE(ZRHS(IMICRO))
  ALLOCATE(ZTHS(IMICRO))
  ALLOCATE(ZRHODREF(IMICRO))
  ALLOCATE(ZZT(IMICRO))
  ALLOCATE(ZPRES(IMICRO))
  ALLOCATE(ZEXNREF(IMICRO))
  ALLOCATE(ZSIGMA_RC(IMICRO))
  ALLOCATE(ZCF(IMICRO))
  DO JL = 1, IMICRO   
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
    IF (KRR == 7) ZRHT(JL) = PRHT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
    IF (HSUBG_AUCV == 'SIGM') THEN
         ZSIGMA_RC(JL) = PSIGS(I1(JL),I2(JL),I3(JL)) * 2.
      ELSE IF (HSUBG_AUCV == 'CLFR') THEN
         ZCF(JL) = PCLDFR(I1(JL),I2(JL),I3(JL))
    END IF
!
    ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
    ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
    ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
    ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
    ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
    ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
    IF (KRR == 7) ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
    ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
  ENDDO
!
  ALLOCATE(ZRHOCOR(IMICRO))
  ZRHOCOR(:) = (ZRHO00 / ZRHODREF(:))**XCEXVT
!
  ALLOCATE(ZZW(IMICRO))
  ALLOCATE(ZLSFACT(IMICRO))
  ALLOCATE(ZLVFACT(IMICRO))
!
  ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
                                  +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
  ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
  ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)
!
  ALLOCATE(ZUSW(IMICRO))
  ALLOCATE(ZSSI(IMICRO))
!
  ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
  ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                    ! Supersaturation over ice
!
  IF (KRR == 7) THEN
    ALLOCATE(ZRSMIN_ELEC(7))
  ELSE
    ALLOCATE(ZRSMIN_ELEC(6))
  END IF
  ZRSMIN_ELEC(:) = XRTMIN_ELEC(:) / PTSTEP
!
  ALLOCATE(ZLBDAR(IMICRO))
  ALLOCATE(ZLBDAS(IMICRO))
  ALLOCATE(ZLBDAG(IMICRO))
  IF (KRR == 7) ALLOCATE(ZLBDAH(IMICRO))
  ALLOCATE(ZRDRYG(IMICRO))
  ALLOCATE(ZRWETG(IMICRO))
  ALLOCATE(ZAI(IMICRO))
  ALLOCATE(ZCJ(IMICRO))
  ALLOCATE(ZKA(IMICRO))
  ALLOCATE(ZDV(IMICRO))
!
  IF (KRR == 7) THEN
    ALLOCATE(ZZW1(IMICRO,7))
  ELSE IF(KRR == 6) THEN
    ALLOCATE(ZZW1(IMICRO,6))
  ENDIF
!
  IF (LBU_ENABLE .OR. LLES_CALL) THEN
    ALLOCATE(ZRHODJ(IMICRO))
    DO JL=1,IMICRO
      ZRHODJ(JL) = PRHODJ(I1(JL),I2(JL),I3(JL))
    END DO
  END IF
!
  ALLOCATE( ZECT(IMICRO) )
  ALLOCATE( ZERT(IMICRO) )
  ALLOCATE( ZEIT(IMICRO) )
  ALLOCATE( ZEST(IMICRO) )
  ALLOCATE( ZEGT(IMICRO) )
  IF ( KRR == 7 ) ALLOCATE(ZEHT(IMICRO))
  ALLOCATE( ZECS(IMICRO) )
  ALLOCATE( ZERS(IMICRO) )
  ALLOCATE( ZEIS(IMICRO) )
  ALLOCATE( ZESS(IMICRO) )
  ALLOCATE( ZEGS(IMICRO) )
  IF ( KRR == 7 ) ALLOCATE(ZEHS(IMICRO))
  ALLOCATE( ZQPIT(IMICRO) )           
  ALLOCATE( ZQNIT(IMICRO) )           
  ALLOCATE( ZQCT(IMICRO) ) 
  ALLOCATE( ZQRT(IMICRO) ) 
  ALLOCATE( ZQIT(IMICRO) ) 
  ALLOCATE( ZQST(IMICRO) ) 
  ALLOCATE( ZQGT(IMICRO) )
  IF ( KRR == 7 ) ALLOCATE(ZQHT(IMICRO))
  ALLOCATE( ZQPIS(IMICRO) )           
  ALLOCATE( ZQNIS(IMICRO) )           
  ALLOCATE( ZQCS(IMICRO) ) 
  ALLOCATE( ZQRS(IMICRO) ) 
  ALLOCATE( ZQIS(IMICRO) )
  ALLOCATE( ZQSS(IMICRO) ) 
  ALLOCATE( ZQGS(IMICRO) ) 
  IF ( KRR == 7 ) ALLOCATE(ZQHS(IMICRO))
!
  IF (CNI_CHARGING == 'GARDI') THEN
    ALLOCATE( ZDELTALWC(IMICRO) )
    ALLOCATE( ZFT(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TAKAH' .OR.                              &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
    ALLOCATE( ZEW(IMICRO) )
  END IF

  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
    ALLOCATE( ZLWCC(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'TEEWC') THEN
    ALLOCATE( ZDQLWC(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC' )  THEN
    ALLOCATE( ZSAUNSK(IMICRO) )
    ALLOCATE( ZSAUNIM(IMICRO) )
    ALLOCATE( ZSAUNIN(IMICRO) )
    ALLOCATE( ZSAUNSM(IMICRO) )
    ALLOCATE( ZSAUNSN(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'SAP98' .OR.                              &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
    ALLOCATE( ZFQIAGGS(IMICRO) )
    ALLOCATE( ZFQIDRYGBS(IMICRO) )
    ALLOCATE( ZLBQSDRYGB1S(IMICRO) )
    ALLOCATE( ZLBQSDRYGB2S(IMICRO) )
    ALLOCATE( ZLBQSDRYGB3S(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAP98' .OR.                              &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ALLOCATE( ZRAR_CRIT(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'TERAR' .OR. &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ALLOCATE( ZVGMEAN(IMICRO) )
    ALLOCATE( ZVSMEAN(IMICRO) )
    ALLOCATE( ZRAR(IMICRO) )
    ALLOCATE( ZDQRAR_IS(IMICRO) )
    ALLOCATE( ZDQRAR_IG(IMICRO) )
    ALLOCATE( ZDQRAR_SG(IMICRO) )
    ALLOCATE( ZSAUNIM_IS(IMICRO) )
    ALLOCATE( ZSAUNIN_IS(IMICRO) )
    ALLOCATE( ZSAUNIM_IG(IMICRO) )
    ALLOCATE( ZSAUNIN_IG(IMICRO) )
    ALLOCATE( ZSAUNSK_SG(IMICRO) )
    ALLOCATE( ZSAUNSM_SG(IMICRO) )
    ALLOCATE( ZSAUNSN_SG(IMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAP98' .OR. &
      CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'GARDI' .OR. CNI_CHARGING == 'BSMP1' .OR. &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TEEWC' .OR. &
      CNI_CHARGING == 'TERAR') THEN
    ALLOCATE( ZAUX1(IMICRO) )
    ALLOCATE( ZLIMIT(IMICRO) )
  END IF
!
  IF (LINDUCTIVE) THEN
    ALLOCATE( ZEFIELDW(IMICRO) )
    ALLOCATE( ZRATE_IND(IMICRO) )
    ALLOCATE( GIND(IMICRO) )
  END IF
!
  IF (LIAGGS_LATHAM) THEN
    ALLOCATE( ZEFIELDU(IMICRO) )
    ALLOCATE( ZEFIELDV(IMICRO) )
    IF (.NOT.ALLOCATED(ZEFIELDW)) ALLOCATE( ZEFIELDW(IMICRO) )
  END IF
  ALLOCATE( ZLATHAMIAGGS(IMICRO) )
!
  ALLOCATE( ZWQ1(IMICRO,10) )
  ALLOCATE( ZWQ3(IMICRO) )
  ALLOCATE( ZWQ4(IMICRO) )
  ALLOCATE( ZCOLIS(IMICRO) )
  ALLOCATE( ZCOLIG(IMICRO) )
  ALLOCATE( ZCOLSG(IMICRO) )
  ALLOCATE( GELEC(IMICRO,4) )
  GELEC(:,:) = .FALSE.
!
  DO JL = 1, IMICRO
    IF (LINDUCTIVE) ZEFIELDW(JL) = XEFIELDW(I1(JL), I2(JL), I3(JL))
    IF (LIAGGS_LATHAM) THEN
      ZEFIELDU(JL) = XEFIELDU(I1(JL), I2(JL), I3(JL))
      ZEFIELDV(JL) = XEFIELDV(I1(JL), I2(JL), I3(JL))
      IF (.NOT.LINDUCTIVE ) ZEFIELDW(JL) = XEFIELDW(I1(JL), I2(JL), I3(JL))
    END IF
!
    ZQPIT(JL) = PQPIT(I1(JL), I2(JL), I3(JL))
    ZQNIT(JL) = PQNIT(I1(JL), I2(JL), I3(JL))
    ZQCT(JL) = PQCT(I1(JL), I2(JL), I3(JL))
    ZQRT(JL) = PQRT(I1(JL), I2(JL), I3(JL))
    ZQIT(JL) = PQIT(I1(JL), I2(JL), I3(JL))
    ZQST(JL) = PQST(I1(JL), I2(JL), I3(JL))
    ZQGT(JL) = PQGT(I1(JL), I2(JL), I3(JL))
    IF (KRR == 7) ZQHT(JL) = PQHT(I1(JL), I2(JL), I3(JL))
!
    ZQPIS(JL) = PQPIS(I1(JL), I2(JL), I3(JL))
    ZQNIS(JL) = PQNIS(I1(JL), I2(JL), I3(JL))
    ZQCS(JL) = PQCS(I1(JL), I2(JL), I3(JL))
    ZQRS(JL) = PQRS(I1(JL), I2(JL), I3(JL))
    ZQIS(JL) = PQIS(I1(JL), I2(JL), I3(JL))
    ZQSS(JL) = PQSS(I1(JL), I2(JL), I3(JL))
    ZQGS(JL) = PQGS(I1(JL), I2(JL), I3(JL))
    IF (KRR == 7) ZQHS(JL) = PQHS(I1(JL), I2(JL), I3(JL))
  ENDDO
!
!
!*       2.3    Update the parameter e in the charge-diameter relation
!
  IF (KRR == 7) THEN
    CALL COMPUTE_LBDA(ZRRT, ZRST, ZRGT, ZRH=ZRHT)
    ZTSPLITR = 1.
    CALL ELEC_UPDATE_QD(ZTSPLITR, ZERT, ZEIT, ZEST, ZEGT, ZQRT, ZQIT, ZQST, ZQGT, &
                        ZRRT, ZRIT, ZRST, ZRGT,                                   &
                        ZEH=ZEHT, ZQH=ZQHT, ZRH=ZRHT, ZEC=ZECT, ZQC=ZQCT, ZRC=ZRCT)
    ZTSPLITR = PTSTEP
    CALL ELEC_UPDATE_QD(ZTSPLITR, ZERS, ZEIS, ZESS, ZEGS, ZQRS, ZQIS, ZQSS, ZQGS, &
                        ZRRS, ZRIS, ZRSS, ZRGS,                                   &
                        ZEH=ZEHS, ZQH=ZQHS, ZRH=ZRHS, ZEC=ZECS, ZQC=ZQCS, ZRC=ZRCS)
  ELSE
    CALL COMPUTE_LBDA(ZRRT, ZRST, ZRGT)
    ZTSPLITR = 1.
    CALL ELEC_UPDATE_QD(ZTSPLITR, ZERT, ZEIT, ZEST, ZEGT, ZQRT, ZQIT, ZQST, ZQGT, &
                        ZRRT, ZRIT, ZRST, ZRGT, ZEC=ZECT, ZQC=ZQCT, ZRC=ZRCT)
    ZTSPLITR = PTSTEP
    CALL ELEC_UPDATE_QD(ZTSPLITR, ZERS, ZEIS, ZESS, ZEGS, ZQRS, ZQIS, ZQSS, ZQGS, &
                        ZRRS, ZRIS, ZRSS, ZRGS, ZEC=ZECS, ZQC=ZQCS, ZRC=ZRCS)
  END IF
!
!
!*   	 2.4 	Initialization for the non-inductive charging process	
!
  CALL ELEC_INI_NI_PROCESS 
!
!
!*	 2.5	Compute the slow cold process sources
!
  CALL RAIN_ICE_ELEC_SLOW
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTES THE SLOW WARM PROCESS SOURCES
!   	        --------------------------------------
!
  IF( OWARM ) THEN    !  Check if the formation of the raindrops by the slow
                      !  warm processes is allowed
    PEVAP3D(:,:,:)= 0.
    CALL RAIN_ICE_ELEC_WARM
  END IF
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
  CALL RAIN_ICE_ELEC_FAST_RS
!
!-------------------------------------------------------------------------------
!
!*       5.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!               ----------------------------------------------
!
  CALL RAIN_ICE_ELEC_FAST_RG
!
!-------------------------------------------------------------------------------
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
  IF ( KRR == 7 ) THEN
    CALL RAIN_ICE_ELEC_FAST_RH
  END IF
!
!-------------------------------------------------------------------------------
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!
  CALL RAIN_ICE_ELEC_FAST_RI
!
!
!-------------------------------------------------------------------------------
!
!*       8.     UPDATE MIXING 3D RATIOS AND VOLUMETRIC CHARGE CONCENTRATIONS
!		------------------------------------------------------------
!
!*       8.1     Update the mixing ratio
!
  DO JL=1,IMICRO
    PRVS(I1(JL),I2(JL),I3(JL)) = ZRVS(JL)
    PRCS(I1(JL),I2(JL),I3(JL)) = ZRCS(JL)
    PRRS(I1(JL),I2(JL),I3(JL)) = ZRRS(JL)
    PRIS(I1(JL),I2(JL),I3(JL)) = ZRIS(JL)
    PRSS(I1(JL),I2(JL),I3(JL)) = ZRSS(JL)
    PRGS(I1(JL),I2(JL),I3(JL)) = ZRGS(JL)
    PTHS(I1(JL),I2(JL),I3(JL)) = ZTHS(JL)
    PCIT(I1(JL),I2(JL),I3(JL)) = ZCIT(JL)
  END DO
  IF ( KRR == 7 ) THEN
    DO JL=1,IMICRO
      PRHS(I1(JL),I2(JL),I3(JL)) = ZRHS(JL)
    END DO
  END IF
!
!
!*	8.2	Compute the volumetric charge concentration
!
  DO JL=1,IMICRO
    PQPIS(I1(JL),I2(JL),I3(JL)) = ZQPIS(JL)
    PQNIS(I1(JL),I2(JL),I3(JL)) = ZQNIS(JL)
    PQCS (I1(JL),I2(JL),I3(JL)) = ZQCS(JL)
    PQRS (I1(JL),I2(JL),I3(JL)) = ZQRS(JL)
    PQIS (I1(JL),I2(JL),I3(JL)) = ZQIS(JL)
    PQSS (I1(JL),I2(JL),I3(JL)) = ZQSS(JL)
    PQGS (I1(JL),I2(JL),I3(JL)) = ZQGS(JL)
  END DO
  IF ( KRR == 7 ) THEN
    DO JL=1,IMICRO
      PQHS(I1(JL),I2(JL),I3(JL)) = ZQHS(JL)
    END DO
  END IF
!
!
!*      8.3     Deallocate
!
  DEALLOCATE(ZZW1)
  DEALLOCATE(ZDV)
  DEALLOCATE(ZCJ)
  DEALLOCATE(ZRDRYG)
  DEALLOCATE(ZRWETG)
  DEALLOCATE(ZLBDAG)
  IF ( KRR == 7 ) DEALLOCATE(ZLBDAH)
  DEALLOCATE(ZLBDAS)
  DEALLOCATE(ZLBDAR)
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZLVFACT)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZRHODREF)
  DEALLOCATE(ZRHOCOR)
  DEALLOCATE(ZZT)
  IF(LBU_ENABLE .OR. LLES_CALL) DEALLOCATE(ZRHODJ)
  DEALLOCATE(ZTHS)
  IF ( KRR == 7 ) DEALLOCATE(ZRHS)
  DEALLOCATE(ZRGS)
  DEALLOCATE(ZRSS)
  DEALLOCATE(ZRIS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZRGT)
  IF ( KRR == 7 ) DEALLOCATE(ZRHT)
  DEALLOCATE(ZRST)
  DEALLOCATE(ZRIT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZAI)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZKA)
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZSIGMA_RC)
  DEALLOCATE(ZCF)
!
  DEALLOCATE( ZECT )
  DEALLOCATE( ZERT )
  DEALLOCATE( ZEIT )
  DEALLOCATE( ZEST )
  DEALLOCATE( ZEGT )
  IF ( KRR == 7 ) DEALLOCATE(ZEHT)
  DEALLOCATE( ZECS )
  DEALLOCATE( ZERS )
  DEALLOCATE( ZEIS )
  DEALLOCATE( ZESS )
  DEALLOCATE( ZEGS )
  IF ( KRR == 7 ) DEALLOCATE(ZEHS)
  DEALLOCATE( ZQPIT ) 
  DEALLOCATE( ZQNIT ) 
  DEALLOCATE( ZQCT ) 
  DEALLOCATE( ZQRT ) 
  DEALLOCATE( ZQIT )
  DEALLOCATE( ZQST ) 
  DEALLOCATE( ZQGT ) 
  IF ( KRR == 7 ) DEALLOCATE(ZQHT) 
  DEALLOCATE( ZQPIS ) 
  DEALLOCATE( ZQNIS ) 
  DEALLOCATE( ZQCS ) 
  DEALLOCATE( ZQRS ) 
  DEALLOCATE( ZQIS )
  DEALLOCATE( ZQSS )  
  DEALLOCATE( ZQGS ) 
  IF ( KRR == 7 ) DEALLOCATE(ZQHS)
  DEALLOCATE( ZWQ1 )
  DEALLOCATE( ZWQ3 )
  DEALLOCATE( ZWQ4 )
  DEALLOCATE( ZCOLIS )
  DEALLOCATE( ZCOLIG )
  DEALLOCATE( ZCOLSG )
  DEALLOCATE( ZRSMIN_ELEC)
  DEALLOCATE( GELEC )
  IF (ALLOCATED( ZDELTALWC )) DEALLOCATE( ZDELTALWC )   
  IF (ALLOCATED( ZLWCC ))     DEALLOCATE( ZLWCC )
  IF (ALLOCATED( ZFT ))       DEALLOCATE( ZFT )
  IF (ALLOCATED( ZEW ))       DEALLOCATE( ZEW )
  IF (ALLOCATED( ZSAUNSK ))   DEALLOCATE( ZSAUNSK )
  IF (ALLOCATED( ZSAUNIM ))   DEALLOCATE( ZSAUNIM )
  IF (ALLOCATED( ZSAUNIN ))   DEALLOCATE( ZSAUNIN )
  IF (ALLOCATED( ZSAUNSM ))   DEALLOCATE( ZSAUNSM )
  IF (ALLOCATED( ZSAUNSN ))   DEALLOCATE( ZSAUNSN )
  IF (ALLOCATED( ZVGMEAN ))   DEALLOCATE( ZVGMEAN )
  IF (ALLOCATED( ZRAR ))      DEALLOCATE( ZRAR )
  IF (ALLOCATED( ZRAR_CRIT )) DEALLOCATE( ZRAR_CRIT )
  IF (ALLOCATED( ZSAUNIM_IS ))   DEALLOCATE( ZSAUNIM_IS )
  IF (ALLOCATED( ZSAUNIN_IS ))   DEALLOCATE( ZSAUNIN_IS )
  IF (ALLOCATED( ZFQIAGGS ))   DEALLOCATE( ZFQIAGGS )
  IF (ALLOCATED( ZFQIDRYGBS ))   DEALLOCATE( ZFQIDRYGBS )
  IF (ALLOCATED( ZLBQSDRYGB1S ))   DEALLOCATE( ZLBQSDRYGB1S )
  IF (ALLOCATED( ZLBQSDRYGB2S ))   DEALLOCATE( ZLBQSDRYGB2S )
  IF (ALLOCATED( ZLBQSDRYGB3S ))   DEALLOCATE( ZLBQSDRYGB3S )
  IF (ALLOCATED( ZSAUNIM_IG ))   DEALLOCATE( ZSAUNIM_IG )
  IF (ALLOCATED( ZSAUNIN_IG ))   DEALLOCATE( ZSAUNIN_IG )
  IF (ALLOCATED( ZSAUNSK_SG ))   DEALLOCATE( ZSAUNSK_SG )
  IF (ALLOCATED( ZSAUNSM_SG ))   DEALLOCATE( ZSAUNSM_SG )
  IF (ALLOCATED( ZSAUNSN_SG ))   DEALLOCATE( ZSAUNSN_SG )
  IF (ALLOCATED( ZDQLWC ))       DEALLOCATE( ZDQLWC )
  IF (ALLOCATED( ZDQRAR_IS ))    DEALLOCATE( ZDQRAR_IS )
  IF (ALLOCATED( ZDQRAR_IG ))    DEALLOCATE( ZDQRAR_IG )
  IF (ALLOCATED( ZDQRAR_SG ))    DEALLOCATE( ZDQRAR_SG )
  IF (ALLOCATED( ZAUX1 ))     DEALLOCATE( ZAUX1 )
  IF (ALLOCATED( ZLIMIT ))    DEALLOCATE( ZLIMIT )
  IF (ALLOCATED( ZEFIELDW ))  DEALLOCATE( ZEFIELDW )
  IF (ALLOCATED( ZRATE_IND )) DEALLOCATE( ZRATE_IND )
  IF (ALLOCATED( GIND ))      DEALLOCATE( GIND )
  IF (ALLOCATED( ZEFIELDU ))  DEALLOCATE( ZEFIELDU )
  IF (ALLOCATED( ZEFIELDV ))  DEALLOCATE( ZEFIELDV )
  DEALLOCATE( ZLATHAMIAGGS )
!
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
!*       8.1    time splitting loop initialization
!
ZTSPLITR = PTSTEP / REAL(KSPLITR)
!
!
IF (CSEDIM == 'STAT') THEN
! not yet developped for electricity !!!
  CALL RAIN_ICE_SEDIMENTATION_STAT
ELSE
  CALL RAIN_ICE_ELEC_SEDIMENTATION_SPLIT
END IF
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE RAIN_ICE_ELEC_SEDIMENTATION_SPLIT
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
INTEGER , DIMENSION(SIZE(GSEDIMC)) :: IC1,IC2,IC3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMR)) :: IR1,IR2,IR3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMI)) :: II1,II2,II3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMS)) :: IS1,IS2,IS3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMG)) :: IG1,IG2,IG3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMH)) :: IH1,IH2,IH3 ! Used to replace the COUNT
INTEGER   :: ILENALLOCC,ILENALLOCR,ILENALLOCI,ILENALLOCS,ILENALLOCG,ILENALLOCH
INTEGER   :: ILISTLENC,ILISTLENR,ILISTLENI,ILISTLENS,ILISTLENG,ILISTLENH
INTEGER, ALLOCATABLE :: ILISTR(:),ILISTC(:),ILISTI(:),ILISTS(:),ILISTG(:),ILISTH(:)
! Optimization for NEC
!INTEGER, SAVE :: IOLDALLOCC = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
!INTEGER, SAVE :: IOLDALLOCR = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
!INTEGER, SAVE :: IOLDALLOCI = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
!INTEGER, SAVE :: IOLDALLOCS = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
!INTEGER, SAVE :: IOLDALLOCG = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
!INTEGER, SAVE :: IOLDALLOCH = SIZE(PEXNREF,1)*SIZE(PEXNREF,2)*SIZE(PEXNREF,3)/10
INTEGER, SAVE :: IOLDALLOCC = 6000
INTEGER, SAVE :: IOLDALLOCR = 6000
INTEGER, SAVE :: IOLDALLOCI = 6000
INTEGER, SAVE :: IOLDALLOCS = 6000
INTEGER, SAVE :: IOLDALLOCG = 6000
INTEGER, SAVE :: IOLDALLOCH = 6000
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D !  droplet condensation
INTEGER, DIMENSION(:), ALLOCATABLE :: ZCIS
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZF0, ZF1, ZCOR
REAL :: ZBEARDCOEFR, ZBEARDCOEFI, ZBEARDCOEFS, ZBEARDCOEFG
REAL :: ZVR, ZVI, ZVS, ZVG, ZETA0, ZK, ZRE0
! For rain, ice, snow and graupel particles, Take into account the 
! effects of altitude and electrical force on terminal fallspeed
! (from Beard, JAS 1980, 37,1363-1374) 
!
!-------------------------------------------------------------------------------
!
!        O. Initialization for sedimentation                  
!
  if ( lbudget_rc .and. osedic ) &
                    call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( osedic ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'SEDI', pqcs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'SEDI', pqrs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'SEDI', pqis(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'SEDI', pqss(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'SEDI', pqgs(:, :, :) * prhodj(:, :, :) )
    if ( krr == 7 ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'SEDI', pqhs(:, :, :) * prhodj(:, :, :) )
  end if

  IF (OSEDIC) PINPRC (:,:) = 0.
  PINPRR (:,:) = 0.
  PINPRR3D (:,:,:) = 0.
  PINPRS (:,:) = 0.
  PINPRG (:,:) = 0.
  IF ( KRR == 7 ) PINPRH (:,:) = 0.
!
  ZT (:,:,:) = ZT (:,:,:) - XTT    !ZT from RAIN_ICE_ELEC_NUCLEATION
  ZETA0 = (1.718 + 0.0049*(XTHVREFZ(IKB) -XTT))
  WHERE (ZT (:,:,:) >= 0.0)
    ZF0(:,:,:) = ZETA0 / (1.718 + 0.0049*ZT(:,:,:))
  ELSEWHERE
    ZF0(:,:,:) = ZETA0 / (1.718 + 0.0049*ZT(:,:,:) - 1.2E-5*ZT(:,:,:)*ZT(:,:,:))
  END WHERE
!
  ZF1(:,:,:) = SQRT(ZRHO00/PRHODREF(:,:,:))
  ZCOR(:,:,:) = (PRHODREF(:,:,:)/ZRHO00)**XCEXVT   ! to eliminate Foote-duToit correction
!
  ZVR = (ZRHO00/ZETA0) * XCR * MOMG(XALPHAR,XNUR,XBR+XDR) / MOMG(XALPHAR,XNUR,XBR)
  ZVI = (ZRHO00/ZETA0) * 2.1E5 * MOMG(XALPHAI,XNUI,3.285) / MOMG(XALPHAI,XNUI,1.7)  ! Columns
  ZVS = (ZRHO00/ZETA0) * XCS * MOMG(XALPHAS,XNUS,XBS+XDS) / MOMG(XALPHAS,XNUS,XBS)
  ZVG = (ZRHO00/ZETA0) * XCG * MOMG(XALPHAG,XNUG,XBG+XDG) / MOMG(XALPHAG,XNUG,XBG)
!
!*       1. Parameters for cloud sedimentation
!
  IF (OSEDIC) THEN
    ZRAY(:,:,:)    = 0.
    ZLBC(:,:,:)    = XLBC(1)
    ZFSEDC(:,:,:)  = XFSEDC(1)
    ZCONC3D(:,:,:) = XCONC_LAND
    ZCONC_TMP(:,:) = XCONC_LAND
    IF (PRESENT(PSEA)) THEN
      ZCONC_TMP(:,:) = PSEA(:,:) * XCONC_SEA + (1. - PSEA(:,:)) * XCONC_LAND
      DO JK = IKB, IKE
        ZLBC(:,:,JK)    = PSEA(:,:) * XLBC(2)   + (1. - PSEA(:,:)) * XLBC(1)
        ZFSEDC(:,:,JK)  = PSEA(:,:) * XFSEDC(2) + (1. - PSEA(:,:)) * XFSEDC(1)
        ZFSEDC(:,:,JK)  = MAX(MIN(XFSEDC(1),XFSEDC(2)),ZFSEDC(:,:,JK))
        ZCONC3D(:,:,JK) = (1. - PTOWN(:,:)) * ZCONC_TMP(:,:) + PTOWN(:,:) * XCONC_URBAN
        ZRAY(:,:,JK)    = 0.5 * ((1. - PSEA(:,:)) * MOMG(XALPHAC, XNUC, 1.0) + &
                                       PSEA(:,:)  * MOMG(XALPHAC2, XNUC2, 1.0) )
      END DO
    ELSE
        ZCONC3D(:,:,:) = XCONC_LAND
        ZRAY(:,:,:)    = 0.5 * MOMG(XALPHAC, XNUC, 1.0)
    END IF
    ZRAY(:,:,:) = MAX(1.,ZRAY(:,:,:))
    ZLBC(:,:,:) = MAX(MIN(XLBC(1),XLBC(2)),ZLBC(:,:,:))
  ENDIF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately

  ZRTMIN(:) = XRTMIN(:) / PTSTEP
  IF (OSEDIC) GSEDIMC(:,:,:) = .FALSE.
  GSEDIMR(:,:,:) = .FALSE.
  GSEDIMI(:,:,:) = .FALSE.
  GSEDIMS(:,:,:) = .FALSE.
  GSEDIMG(:,:,:) = .FALSE.
  IF (KRR == 7) GSEDIMH(:,:,:) = .FALSE.
!
  IF (OSEDIC) ILENALLOCC = 0
  ILENALLOCR = 0
  ILENALLOCI = 0
  ILENALLOCS = 0
  ILENALLOCG = 0
  IF ( KRR == 7 ) ILENALLOCH = 0
!
! ZPiS = Specie i source creating during the current time step
! PRiS = Source of the previous time step
!
  IF (OSEDIC) THEN
    ZPRCS(:,:,:) = 0.0
    ZPRCS(:,:,:) = PRCS(:,:,:) - PRCT(:,:,:) / PTSTEP
    PRCS(:,:,:)  = PRCT(:,:,:) / PTSTEP
  END IF
  ZPRRS(:,:,:) = 0.0
  ZPRSS(:,:,:) = 0.0
  ZPRGS(:,:,:) = 0.0
  IF (KRR == 7) ZPRHS(:,:,:) = 0.0
!
  ZPRRS(:,:,:) = PRRS(:,:,:) - PRRT(:,:,:) / PTSTEP
  ZPRSS(:,:,:) = PRSS(:,:,:) - PRST(:,:,:) / PTSTEP
  ZPRGS(:,:,:) = PRGS(:,:,:) - PRGT(:,:,:) / PTSTEP
  IF (KRR == 7) ZPRHS(:,:,:) = PRHS(:,:,:) - PRHT(:,:,:) / PTSTEP
  PRRS(:,:,:)  = PRRT(:,:,:) / PTSTEP
  PRSS(:,:,:)  = PRST(:,:,:) / PTSTEP
  PRGS(:,:,:)  = PRGT(:,:,:) / PTSTEP
  IF (KRR == 7) PRHS(:,:,:) = PRHT(:,:,:) / PTSTEP
  ZPQRS(:,:,:) = 0.0
  ZPQSS(:,:,:) = 0.0
  ZPQGS(:,:,:) = 0.0
  IF (KRR == 7) ZPQHS(:,:,:) = 0.0
!
  ZPQRS(:,:,:) = PQRS(:,:,:) - PQRT(:,:,:) / PTSTEP
  ZPQSS(:,:,:) = PQSS(:,:,:) - PQST(:,:,:) / PTSTEP
  ZPQGS(:,:,:) = PQGS(:,:,:) - PQGT(:,:,:) / PTSTEP
  IF (KRR == 7) ZPQHS(:,:,:) = PQHS(:,:,:) - PQHT(:,:,:) / PTSTEP
  PQRS(:,:,:)  = PQRT(:,:,:) / PTSTEP
  PQSS(:,:,:)  = PQST(:,:,:) / PTSTEP
  PQGS(:,:,:)  = PQGT(:,:,:) / PTSTEP
  IF (KRR == 7) PQHS(:,:,:) = PQHT(:,:,:) / PTSTEP
!
! PRiS = Source of the previous time step + source created during the subtime
! step
!
  DO JN = 1, KSPLITR
    IF(JN == 1) THEN
      IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:) / KSPLITR
      PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:) / KSPLITR
      PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:) / KSPLITR
      PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:) / KSPLITR
      IF (KRR == 7) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:) / KSPLITR
      PQRS(:,:,:) = PQRS(:,:,:) + ZPQRS(:,:,:) / KSPLITR
      PQSS(:,:,:) = PQSS(:,:,:) + ZPQSS(:,:,:) / KSPLITR
      PQGS(:,:,:) = PQGS(:,:,:) + ZPQGS(:,:,:) / KSPLITR
      IF (KRR == 7) PQHS(:,:,:) = PQHS(:,:,:) + ZPQHS(:,:,:) / KSPLITR
      DO JK = IKB, IKE
        ZW(:,:,JK) = ZTSPLITR / (PRHODREF(:,:,JK) * (PZZ(:,:,JK+1) - PZZ(:,:,JK)))
      END DO
    ELSE
      IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:) * ZTSPLITR
      PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:) * ZTSPLITR
      PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:) * ZTSPLITR
      PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:) * ZTSPLITR
      IF (KRR == 7) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:) * ZTSPLITR
      PQRS(:,:,:) = PQRS(:,:,:) + ZPQRS(:,:,:) * ZTSPLITR
      PQSS(:,:,:) = PQSS(:,:,:) + ZPQSS(:,:,:) * ZTSPLITR
      PQGS(:,:,:) = PQGS(:,:,:) + ZPQGS(:,:,:) * ZTSPLITR
      IF (KRR == 7) PQHS(:,:,:) = PQHS(:,:,:) + ZPQHS(:,:,:) * ZTSPLITR
    END IF
 !
    IF (OSEDIC) GSEDIMC(IIB:IIE,IJB:IJE,IKB:IKE) =                &
                     PRCS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(2)
    GSEDIMR(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                     PRRS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(3)
    GSEDIMI(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                     PRIS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(4)
    GSEDIMS(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                     PRSS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(5)
    GSEDIMG(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
                     PRGS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(6)
    IF (KRR == 7) GSEDIMH(IIB:IIE,IJB:IJE,IKB:IKE) =              &
                     PRHS(IIB:IIE,IJB:IJE,IKB:IKE) > ZRTMIN(7)
!
    IF (OSEDIC) ISEDIMC = COUNTJV( GSEDIMC(:,:,:),IC1(:),IC2(:),IC3(:))
    ISEDIMR = COUNTJV( GSEDIMR(:,:,:),IR1(:),IR2(:),IR3(:))
    ISEDIMI = COUNTJV( GSEDIMI(:,:,:),II1(:),II2(:),II3(:))
    ISEDIMS = COUNTJV( GSEDIMS(:,:,:),IS1(:),IS2(:),IS3(:))
    ISEDIMG = COUNTJV( GSEDIMG(:,:,:),IG1(:),IG2(:),IG3(:))
    IF (KRR == 7) ISEDIMH = COUNTJV( GSEDIMH(:,:,:),IH1(:),IH2(:),IH3(:))
!
!*       2.1   for cloud
!
    IF (OSEDIC) THEN
      ZWSED(:,:,:) = 0.
      IF( JN==1 ) PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP
      IF(ISEDIMC >= 1) THEN
        IF (ISEDIMC .GT. ILENALLOCC) THEN
          IF (ILENALLOCC .GT. 0) THEN
            DEALLOCATE (ZRCS, ZRHODREFC, ILISTC, ZWLBDC, ZCONC, ZRCT,  &
                        ZZT, ZPRES, ZRAY1D, ZFSEDC1D, ZWLBDA, ZCC )
          END IF
          ILENALLOCC = MAX (IOLDALLOCC, 2*ISEDIMC )
          IOLDALLOCC = ILENALLOCC
          ALLOCATE(ZRCS(ILENALLOCC), ZRHODREFC(ILENALLOCC), ILISTC(ILENALLOCC), &
          ZWLBDC(ILENALLOCC), ZCONC(ILENALLOCC), ZRCT(ILENALLOCC), ZZT(ILENALLOCC), &
          ZPRES(ILENALLOCC), ZRAY1D(ILENALLOCC), ZFSEDC1D(ILENALLOCC), &
          ZWLBDA(ILENALLOCC), ZCC(ILENALLOCC))
        END IF
!
        DO JL = 1, ISEDIMC
          ZRCS(JL)      = PRCS(IC1(JL),IC2(JL),IC3(JL))
          ZRHODREFC(JL) = PRHODREF(IC1(JL),IC2(JL),IC3(JL))
          ZWLBDC(JL)    = ZLBC(IC1(JL),IC2(JL),IC3(JL))
          ZCONC(JL)     = ZCONC3D(IC1(JL),IC2(JL),IC3(JL))
          ZRCT(JL)      = PRCT(IC1(JL),IC2(JL),IC3(JL))
          ZZT(JL)       = PTHT(IC1(JL),IC2(JL),IC3(JL))
          ZPRES(JL)     = PPABST(IC1(JL),IC2(JL),IC3(JL))
          ZRAY1D(JL)    = ZRAY(IC1(JL),IC2(JL),IC3(JL))
          ZFSEDC1D(JL)  = ZFSEDC(IC1(JL),IC2(JL),IC3(JL))
        END DO
!
        ILISTLENC = 0
        DO JL = 1, ISEDIMC
          IF(ZRCS(JL) .GT. ZRTMIN(2)) THEN
            ILISTLENC = ILISTLENC + 1
            ILISTC(ILISTLENC) = JL
          END IF
        END DO
        DO JJ = 1, ILISTLENC
          JL = ILISTC(JJ)
          IF (ZRCS(JL) .GT. ZRTMIN(2) .AND. ZRCT(JL) .GT. XRTMIN(2)) THEN
            ZWLBDC(JL) = ZWLBDC(JL) * ZCONC(JL) / (ZRHODREFC(JL) * ZRCT(JL))
            ZWLBDC(JL) = ZWLBDC(JL)**XLBEXC 
            ZRAY1D(JL) = ZRAY1D(JL) / ZWLBDC(JL) !! ZRAY : mean diameter=M(1)/2
            ZZT(JL)    = ZZT(JL) * (ZPRES(JL) / XP00)**(XRD/XCPD)
            ZWLBDA(JL) = 6.6E-8 * (101325. / ZPRES(JL)) * (ZZT(JL) / 293.15)
            ZCC(JL)    = XCC * (1. + 1.26 * ZWLBDA(JL) / ZRAY1D(JL)) !! XCC modified for cloud
            ZWSED (IC1(JL),IC2(JL),IC3(JL)) = ZRHODREFC(JL)**(-XCEXVT +1 ) *   &
                  ZWLBDC(JL)**(-XDC) * ZCC(JL) * ZFSEDC1D(JL) * ZRCS(JL)
          END IF
        END DO
      END IF
      DO JK = IKB, IKE
        PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
      END DO
      PINPRC(:,:) = PINPRC(:,:) + ZWSED(:,:,IKB) / XRHOLW / KSPLITR 
      IF(JN == KSPLITR) THEN
        PRCS(:,:,:) = PRCS(:,:,:) / PTSTEP
      END IF
    END IF
!
!*       2.2   for rain
!
    IF( JN==1 ) PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP
    IF (JN == 1) PQRS(:,:,:) = PQRS(:,:,:) * PTSTEP
    ZWSED(:,:,:)  = 0.
    ZWSEDQ(:,:,:) = 0.
    IF( ISEDIMR >= 1 ) THEN
      IF ( ISEDIMR .GT. ILENALLOCR ) THEN
        IF ( ILENALLOCR .GT. 0 ) THEN
          DEALLOCATE (ZRRS, ZRHODREFR, ILISTR)
          DEALLOCATE (ZQRS, ZLBDAR, ZERS)
        END IF
        ILENALLOCR = MAX (IOLDALLOCR, 2*ISEDIMR )
        IOLDALLOCR = ILENALLOCR
        ALLOCATE(ZRRS(ILENALLOCR), ZRHODREFR(ILENALLOCR), ILISTR(ILENALLOCR))
        ALLOCATE(ZQRS(ILENALLOCR), ZLBDAR(ILENALLOCR), ZERS(ILENALLOCR))
      END IF
      ZERS(:) = 0.
!
      DO JL = 1, ISEDIMR
        ZRRS(JL)      = PRRS(IR1(JL),IR2(JL),IR3(JL))
        ZRHODREFR(JL) = PRHODREF(IR1(JL),IR2(JL),IR3(JL))
        ZQRS(JL)       = PQRS(IR1(JL),IR2(JL),IR3(JL))
! compute lambda_r and e_r
        IF (ZRRS(JL) > 0.) THEN
          ZLBDAR(JL) = XLBR * (ZRHODREFR(JL) * MAX(ZRRS(JL), ZRTMIN(3)))**XLBEXR
        END IF
        IF (ZRRS(JL) > ZRTMIN(3) .AND. ZLBDAR(JL) > 0.) THEN
          ZERS(JL) = ZRHODREFR(JL) * ZQRS(JL) / (XFQUPDR * ZLBDAR(JL)**(XCXR - XFR))
          ZERS(JL) = SIGN( MIN(ABS(ZERS(JL)), XERMAX), ZERS(JL))
        END IF
      END DO
!
      ILISTLENR = 0
      DO JL = 1, ISEDIMR
        IF(ZRRS(JL) .GT. ZRTMIN(3)) THEN
          ILISTLENR = ILISTLENR + 1
          ILISTR(ILISTLENR) = JL
        END IF
      END DO
      DO JJ = 1, ILISTLENR
        JL = ILISTR(JJ)
        IF (ZRRS(JL) > 0. .AND. LSEDIM_BEARD) THEN
          ZK = 1. - ZQRS(JL) * XEFIELDW(IR1(JL),IR2(JL),IR3(JL)) / (ZRRS(JL)*XG)
          IF (ZK <= 0.0) THEN
            ZBEARDCOEFR = 0.
          ELSE
            ZRE0 = ZVR / ZLBDAR(JL)**(1.+XDR)
            IF (ZRE0 <= 0.2) THEN
              ZBEARDCOEFR = ZF0(IR1(JL),IR2(JL),IR3(JL)) * ZK
            ELSE IF (ZRE0 >= 1000.) THEN
              ZBEARDCOEFR = ZF1(IR1(JL),IR2(JL),IR3(JL)) * SQRT(ZK)
            ELSE
              ZBEARDCOEFR = ZF0(IR1(JL),IR2(JL),IR3(JL)) * ZK +         &
                            (ZF1(IR1(JL),IR2(JL),IR3(JL)) *             &
                            SQRT(ZK)-ZF0(IR1(JL),IR2(JL),IR3(JL))*ZK) * &
                            (1.61+LOG(ZRE0)) / 8.52
            END IF
            ZBEARDCOEFR = ZBEARDCOEFR * ZCOR(IR1(JL),IR2(JL),IR3(JL))
          END IF
        ELSE
          ZBEARDCOEFR = 1.0 ! No "Beard" effect
        END IF
!
        ZWSED(IR1(JL),IR2(JL),IR3(JL)) = ZBEARDCOEFR *                   &
                                         XFSEDR  * ZRRS(JL)**XEXSEDR *   &
                                         ZRHODREFR(JL)**(XEXSEDR-XCEXVT)
!
        IF (ZRRS(JL) > ZRTMIN(3) .AND. ABS(ZERS(JL)) > XERMIN) THEN
          ZWSEDQ(IR1(JL),IR2(JL),IR3(JL)) = ZBEARDCOEFR *                &
                                            XFQSEDR * ZERS(JL) *         &
                                            ZRRS(JL)**XEXQSEDR *         &
                                            ZRHODREFR(JL)**(XEXQSEDR-XCEXVT)
        END IF
      END DO
    END IF
    DO JK = IKB , IKE
      PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1)  - ZWSED(:,:,JK))
      PQRS(:,:,JK) = PQRS(:,:,JK) + ZW(:,:,JK) * (ZWSEDQ(:,:,JK+1) - ZWSEDQ(:,:,JK))
    END DO
    PINPRR(:,:)     = PINPRR(:,:)     + ZWSED(:,:,IKB)  / XRHOLW / KSPLITR
    PINPRR3D(:,:,:) = PINPRR3D(:,:,:) + ZWSED(:,:,:)    / XRHOLW / KSPLITR 
    IF (JN == KSPLITR) THEN
      PRRS(:,:,:) = PRRS(:,:,:) / PTSTEP
      PQRS(:,:,:) = PQRS(:,:,:) / PTSTEP
    END IF
!
!
!*       2.3   for pristine ice
!
    IF (JN == 1) PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
    IF (JN == 1) PQIS(:,:,:) = PQIS(:,:,:) * PTSTEP
    ZWSED(:,:,:)  = 0.
    ZWSEDQ(:,:,:) = 0.
    IF( ISEDIMI >= 1 ) THEN
      IF ( ISEDIMI .GT. ILENALLOCI ) THEN
        IF ( ILENALLOCI .GT. 0 ) THEN
          DEALLOCATE (ZRIS, ZRHODREFI, ILISTI)
          DEALLOCATE (ZQIS, ZEIS, ZCIT, ZCIS, ZLBDAI)
        END IF
        ILENALLOCI = MAX (IOLDALLOCI, 2*ISEDIMI )
        IOLDALLOCI = ILENALLOCI
        ALLOCATE(ZRIS(ILENALLOCI), ZRHODREFI(ILENALLOCI), ILISTI(ILENALLOCI))
        ALLOCATE(ZQIS(ILENALLOCI),       &
                 ZEIS(ILENALLOCI),       &
                 ZCIT(ILENALLOCI),       &
                 ZCIS(ILENALLOCI),       &
                 ZLBDAI(ILENALLOCI))
      END IF
!
      DO JL = 1, ISEDIMI
        ZRIS(JL)      = PRIS(II1(JL),II2(JL),II3(JL))
        ZRHODREFI(JL) = PRHODREF(II1(JL),II2(JL),II3(JL))
        ZQIS(JL)       = PQIS(II1(JL),II2(JL),II3(JL))
        ZCIT(JL)       = PCIT(II1(JL),II2(JL),II3(JL))
        ZEIS(JL) = 0.
! compute e_i
        IF (ZRIS(JL) > ZRTMIN(4) .AND. ZCIT(JL) > 0.0) THEN
          ZEIS(JL) = ZRHODREFI(JL) * ZQIS(JL) / ((ZCIT(JL)**(1 - XEXFQUPDI)) * &
                     XFQUPDI * (ZRHODREFI(JL) * ZRIS(JL))**XEXFQUPDI)
          ZEIS(JL) = SIGN( MIN(ABS(ZEIS(JL)), XEIMAX), ZEIS(JL))
          ZCIS(JL) = XFCI * ZRHODREFI(JL) * ZRIS(JL) * &
                     MAX(0.05E6,                       &
                        -0.15319E6 - 0.021454E6 * ALOG(ZRHODREFI(JL) * ZRIS(JL)))**3
          ZLBDAI(JL) = (2.14E-3 * MOMG(XALPHAI,XNUI,1.7) *     &
                        ZCIS(JL) / (ZRHODREFI(JL) * ZRIS(JL)))**0.588235
        END IF
      END DO
!
      ILISTLENI = 0
      DO JL = 1, ISEDIMI
        IF (ZRIS(JL) .GT. MAX(ZRTMIN(4),1.0E-7 )) THEN ! limitation of the McF&H formula
          ILISTLENI = ILISTLENI + 1
          ILISTI(ILISTLENI) = JL
        END IF
      END DO
      DO JJ = 1, ILISTLENI
        JL = ILISTI(JJ)
        IF (ZRIS(JL) > ZRTMIN(4) .AND. ZCIT(JL) > 0.0 .AND. LSEDIM_BEARD) THEN
          ZK = 1. - ZQIS(JL) * XEFIELDW(II1(JL),II2(JL),II3(JL)) / (ZRIS(JL)*XG)
          IF (ZK <= 0.0) THEN
            ZBEARDCOEFI = 0.
          ELSE
            ZRE0 = ZVI / ZLBDAI(JL)**2.585
            IF (ZRE0 <= 0.2) THEN
              ZBEARDCOEFI = ZF0(II1(JL),II2(JL),II3(JL)) * ZK
            ELSE IF (ZRE0 >= 1000.) THEN
              ZBEARDCOEFI = ZF1(II1(JL),II2(JL),II3(JL)) * SQRT(ZK)
            ELSE
              ZBEARDCOEFI = ZF0(II1(JL),II2(JL),II3(JL)) * ZK +             &
                            (ZF1(II1(JL),II2(JL),II3(JL)) *                 &
                            SQRT(ZK) - ZF0(II1(JL),II2(JL),II3(JL)) * ZK) * &
                            (1.61 + LOG(ZRE0)) / 8.52
            END IF
            ZBEARDCOEFI = ZBEARDCOEFI * ZCOR(II1(JL),II2(JL),II3(JL))
          END IF
        ELSE
          ZBEARDCOEFI = 1.0 ! No "Beard" effect
        END IF
!
        ZWSED(II1(JL),II2(JL),II3(JL))= ZBEARDCOEFI *      &
                             XFSEDI * ZRIS(JL) *           &
                             ZRHODREFI(JL)**(1.0-XCEXVT) * & !    McF&H
                             MAX( 0.05E6,-0.15319E6-0.021454E6* &
                             ALOG(ZRHODREFI(JL)*ZRIS(JL)) )**XEXCSEDI
        IF (ZRIS(JL) .GT. MAX(ZRTMIN(4),1.0E-7) .AND. ABS(ZEIS(JL)) .GT. XEIMIN .AND. &
            ZCIT(JL) .GT. 0. ) THEN
          ZWSEDQ(II1(JL),II2(JL),II3(JL)) = ZBEARDCOEFI *                      &
                     ZCIS(JL)**(1 - XEXQSEDI) * XFQSEDI *                      &
                     ZRIS(JL)**XEXQSEDI * ZRHODREFI(JL)**(XEXQSEDI - XCEXVT) * &
                     ZEIS(JL) * (ZCIT(JL) / ZCIS(JL))**(1.-XFI/XBI)
        END IF
      END DO
    END IF
    DO JK = IKB, IKE
      PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1)  - ZWSED(:,:,JK))
      PQIS(:,:,JK) = PQIS(:,:,JK) + ZW(:,:,JK) * (ZWSEDQ(:,:,JK+1) - ZWSEDQ(:,:,JK))
    END DO
    IF (JN == KSPLITR) THEN
      PRIS(:,:,:) = PRIS(:,:,:) / PTSTEP
      PQIS(:,:,:) = PQIS(:,:,:) / PTSTEP
    END IF
!
!
!*       2.4   for aggregates/snow
!
    IF( JN==1 ) PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
    IF (JN == 1) PQSS(:,:,:) = PQSS(:,:,:) * PTSTEP
    ZWSED(:,:,:)  = 0.
    ZWSEDQ(:,:,:) = 0.
    IF( ISEDIMS >= 1 ) THEN
      IF ( ISEDIMS .GT. ILENALLOCS ) THEN
        IF ( ILENALLOCS .GT. 0 ) THEN
          DEALLOCATE (ZRSS, ZRHODREFS, ILISTS)
          DEALLOCATE (ZQSS, ZESS, ZLBDAS)
        END IF
        ILENALLOCS = MAX(IOLDALLOCS, 2*ISEDIMS )
        IOLDALLOCS = ILENALLOCS
        ALLOCATE(ZRSS(ILENALLOCS), ZRHODREFS(ILENALLOCS), ILISTS(ILENALLOCS))
        ALLOCATE(ZQSS(ILENALLOCS), ZESS(ILENALLOCS), ZLBDAS(ILENALLOCS))
      END IF
!
      DO JL = 1, ISEDIMS
        ZRSS(JL)      = PRSS(IS1(JL),IS2(JL),IS3(JL))
        ZRHODREFS(JL) = PRHODREF(IS1(JL),IS2(JL),IS3(JL))
        ZQSS(JL)       = PQSS(IS1(JL),IS2(JL),IS3(JL))
        ZESS(JL) = 0.
! compute lambda_s and e_s
        IF (ZRSS(JL) > 0.) THEN
          ZLBDAS(JL) = MIN(XLBDAS_MAX, &
                           XLBS * (ZRHODREFS(JL) * MAX(ZRSS(JL), ZRTMIN(5)))**XLBEXS)
        END IF
        IF (ZRSS(JL) > ZRTMIN(5) .AND. ZLBDAS(JL) > 0.) THEN
          ZESS(JL) = ZRHODREFS(JL) * ZQSS(JL) / (XFQUPDS * ZLBDAS(JL)**(XCXS - XFS))
          ZESS(JL) = SIGN( MIN(ABS(ZESS(JL)), XESMAX), ZESS(JL))
        END IF
      END DO
!
      ILISTLENS = 0
      DO JL = 1, ISEDIMS
        IF (ZRSS(JL) .GT. ZRTMIN(5)) THEN
          ILISTLENS = ILISTLENS + 1
          ILISTS(ILISTLENS) = JL
        END IF
      END DO
      DO JJ = 1, ILISTLENS
        JL = ILISTS(JJ)
        IF (ZRSS(JL) > 0. .AND. LSEDIM_BEARD) THEN
          ZK = 1. - ZQSS(JL) * XEFIELDW(IS1(JL),IS2(JL),IS3(JL)) / (ZRSS(JL)*XG)
          IF (ZK <= 0.0) THEN
            ZBEARDCOEFS = 0.
          ELSE
            ZRE0 = ZVS / ZLBDAS(JL)**(1.+XDS)
            IF (ZRE0 <= 0.2) THEN
              ZBEARDCOEFS = ZF0(IS1(JL),IS2(JL),IS3(JL)) * ZK
            ELSE IF (ZRE0 >= 1000.) THEN
              ZBEARDCOEFS = ZF1(IS1(JL),IS2(JL),IS3(JL)) * SQRT(ZK)
            ELSE
              ZBEARDCOEFS = ZF0(IS1(JL),IS2(JL),IS3(JL)) * ZK +            &
                            (ZF1(IS1(JL),IS2(JL),IS3(JL)) *                &
                            SQRT(ZK) -ZF0(IS1(JL),IS2(JL),IS3(JL)) * ZK) * &
                            (1.61 + LOG(ZRE0)) / 8.52
            END IF
            ZBEARDCOEFS = ZBEARDCOEFS * ZCOR(IS1(JL),IS2(JL),IS3(JL))
          END IF
        ELSE
          ZBEARDCOEFS = 1.0 ! No "Beard" effect
        END IF
!
        ZWSED (IS1(JL),IS2(JL),IS3(JL)) = ZBEARDCOEFS *             &
                                      XFSEDS * ZRSS(JL)**XEXSEDS *  &
                                      ZRHODREFS(JL)**(XEXSEDS-XCEXVT)
        IF (ZRSS(JL) .GT. ZRTMIN(5) .AND. ABS(ZESS(JL)) > XESMIN) THEN
          ZWSEDQ(IS1(JL),IS2(JL),IS3(JL)) = ZBEARDCOEFS *         &
                                            XFQSEDS * ZESS(JL) *  &
                                            ZRSS(JL)**XEXQSEDS *  &
                                            ZRHODREFS(JL)**(XEXQSEDS - XCEXVT)
        END IF
      END DO
    END IF
    DO JK = IKB, IKE
      PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1)  - ZWSED(:,:,JK))
      PQSS(:,:,JK) = PQSS(:,:,JK) + ZW(:,:,JK) * (ZWSEDQ(:,:,JK+1) - ZWSEDQ(:,:,JK))
    END DO
    PINPRS(:,:) = PINPRS(:,:) + ZWSED(:,:,IKB)  / XRHOLW / KSPLITR
    IF (JN == KSPLITR) THEN
      PRSS(:,:,:) = PRSS(:,:,:) / PTSTEP
      PQSS(:,:,:) = PQSS(:,:,:) / PTSTEP
    END IF
!
!
!*       2.5   for graupeln
!
    ZWSED(:,:,:)  = 0.
    ZWSEDQ(:,:,:) = 0.
    IF( JN==1 ) PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
    IF (JN == 1) PQGS(:,:,:) = PQGS(:,:,:) * PTSTEP
    IF( ISEDIMG >= 1 ) THEN
      IF ( ISEDIMG .GT. ILENALLOCG ) THEN
        IF ( ILENALLOCG .GT. 0 ) THEN
          DEALLOCATE (ZRGS, ZRHODREFG, ILISTG)
          DEALLOCATE (ZQGS, ZEGS, ZLBDAG)
        END IF
        ILENALLOCG = MAX (IOLDALLOCG, 2*ISEDIMG )
        IOLDALLOCG = ILENALLOCG
        ALLOCATE(ZRGS(ILENALLOCG), ZRHODREFG(ILENALLOCG), ILISTG(ILENALLOCG))
        ALLOCATE(ZQGS(ILENALLOCG), ZEGS(ILENALLOCG), ZLBDAG(ILENALLOCG))
      END IF
!
      DO JL = 1, ISEDIMG
        ZRGS(JL)      = PRGS(IG1(JL),IG2(JL),IG3(JL))
        ZRHODREFG(JL) = PRHODREF(IG1(JL),IG2(JL),IG3(JL))
        ZQGS(JL)      = PQGS(IG1(JL),IG2(JL),IG3(JL))
        ZEGS(JL) = 0.
! compute lambda_g and e_g
        IF (ZRGS(JL) > 0.) THEN
          ZLBDAG(JL) = XLBG * (ZRHODREFG(JL) * MAX(ZRGS(JL), ZRTMIN(6)))**XLBEXG
        END IF
        IF (ZRGS(JL) > ZRTMIN(6) .AND. ZLBDAG(JL) > 0.) THEN
          ZEGS(JL) = ZRHODREFG(JL) * ZQGS(JL) / (XFQUPDG * ZLBDAG(JL)**(XCXG - XFG))
          ZEGS(JL) = SIGN( MIN(ABS(ZEGS(JL)), XEGMAX), ZEGS(JL))
        END IF
      END DO
!
      ILISTLENG = 0
      DO JL = 1, ISEDIMG
        IF (ZRGS(JL) .GT. ZRTMIN(6)) THEN
          ILISTLENG = ILISTLENG + 1
          ILISTG(ILISTLENG) = JL
        END IF
      END DO
      DO JJ = 1, ILISTLENG
        JL = ILISTG(JJ)
        IF (ZRGS(JL) > 0. .AND. LSEDIM_BEARD) THEN
          ZK = 1. - ZQGS(JL) * XEFIELDW(IG1(JL),IG2(JL),IG3(JL)) / (ZRGS(JL)*XG)
          IF (ZK <= 0.0) THEN
            ZBEARDCOEFG = 0.
          ELSE
            ZRE0 = ZVG / ZLBDAG(JL)**(1.+XDG)
            IF (ZRE0 <= 0.2) THEN
              ZBEARDCOEFG = ZF0(IG1(JL),IG2(JL),IG3(JL)) * ZK
            ELSE IF (ZRE0 >= 1000.) THEN
              ZBEARDCOEFG = ZF1(IG1(JL),IG2(JL),IG3(JL)) * SQRT(ZK)
            ELSE
              ZBEARDCOEFG = ZF0(IG1(JL),IG2(JL),IG3(JL)) * ZK +             &
                            (ZF1(IG1(JL),IG2(JL),IG3(JL)) *                 &
                            SQRT(ZK) - ZF0(IG1(JL),IG2(JL),IG3(JL)) * ZK) * &
                            (1.61 + LOG(ZRE0)) / 8.52
            END IF
            ZBEARDCOEFG = ZBEARDCOEFG * ZCOR(IG1(JL),IG2(JL),IG3(JL))
          END IF
        ELSE
          ZBEARDCOEFG = 1.0 ! No "Beard" effect
        END IF
!
        ZWSED (IG1(JL),IG2(JL),IG3(JL))= ZBEARDCOEFG *               &
                                      XFSEDG * ZRGS(JL)**XEXSEDG *   &
                                      ZRHODREFG(JL)**(XEXSEDG-XCEXVT)
        IF (ZRGS(JL) .GT. ZRTMIN(6) .AND. ABS(ZEGS(JL)) > XEGMIN) THEN
          ZWSEDQ(IG1(JL),IG2(JL),IG3(JL)) = ZBEARDCOEFG *         &
                                            XFQSEDG * ZEGS(JL) *  &
                                            ZRGS(JL)**XEXQSEDG *  &
                                            ZRHODREFG(JL)**(XEXQSEDG - XCEXVT)
        END IF
      END DO
    END IF
    DO JK = IKB, IKE
      PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1)  - ZWSED(:,:,JK))
      PQGS(:,:,JK) = PQGS(:,:,JK) + ZW(:,:,JK) * (ZWSEDQ(:,:,JK+1) - ZWSEDQ(:,:,JK))
    END DO
    PINPRG(:,:) = PINPRG(:,:) + ZWSED(:,:,IKB)  / XRHOLW / KSPLITR                             
    IF (JN == KSPLITR) THEN
      PRGS(:,:,:) = PRGS(:,:,:) / PTSTEP
      PQGS(:,:,:) = PQGS(:,:,:) / PTSTEP
    END IF
!
!
!*       2.6   for hail
!
    IF ( KRR == 7 ) THEN
      IF( JN==1 ) PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
      IF (JN == 1) PQHS(:,:,:) = PQHS(:,:,:) * PTSTEP
      ZWSED(:,:,:)  = 0.
      ZWSEDQ(:,:,:) = 0.
      IF( ISEDIMH >= 1 ) THEN
        IF ( ISEDIMH .GT. ILENALLOCH ) THEN
          IF ( ILENALLOCH .GT. 0 ) THEN
            DEALLOCATE (ZRHS, ZRHODREFH, ILISTH)
            DEALLOCATE (ZQHS, ZEHS, ZLBDAH)
          END IF
          ILENALLOCH = MAX(IOLDALLOCH, 2*ISEDIMH )
          IOLDALLOCH = ILENALLOCH
          ALLOCATE(ZRHS(ILENALLOCH), ZRHODREFH(ILENALLOCH), ILISTH(ILENALLOCH))
          ALLOCATE(ZQHS(ILENALLOCH), ZLBDAH(ILENALLOCH), ZEHS(ILENALLOCH))
        END IF
!
        DO JL = 1, ISEDIMH
          ZRHS(JL)      = PRHS(IH1(JL),IH2(JL),IH3(JL))
          ZRHODREFH(JL) = PRHODREF(IH1(JL),IH2(JL),IH3(JL))
          ZQHS(JL)      = PQHS(IH1(JL),IH2(JL),IH3(JL))
          ZEHS(JL) = 0.
! compute lambda_h and e_h
          IF (ZRHS(JL) > 0.) THEN
            ZLBDAH(JL) = XLBH * (ZRHODREFH(JL) * MAX(ZRHS(JL), ZRTMIN(7)))**XLBEXH
          END IF
          IF (ZRHS(JL) > ZRTMIN(7) .AND. ZLBDAH(JL) > 0.) THEN
            ZEHS(JL) = ZRHODREFH(JL) * ZQHS(JL) / (XFQUPDH * ZLBDAH(JL)**(XCXH - XFH))
            ZEHS(JL) = SIGN( MIN(ABS(ZEHS(JL)), XEHMAX), ZEHS(JL))
          END IF
        END DO
!
        ILISTLENH = 0
        DO JL = 1, ISEDIMH
          IF (ZRHS(JL) .GT. ZRTMIN(7)) THEN
            ILISTLENH = ILISTLENH + 1
            ILISTH(ILISTLENH) = JL
          END IF
        END DO
        DO JJ = 1, ILISTLENH
          JL = ILISTH(JJ)
          ZWSED (IH1(JL),IH2(JL),IH3(JL)) = XFSEDH  * ZRHS(JL)**XEXSEDH *   &
                                      ZRHODREFH(JL)**(XEXSEDH-XCEXVT)
          IF (ZRHS(JL) .GT. ZRTMIN(7) .AND. ABS(ZEHS(JL)) > XEHMIN) THEN
            ZWSEDQ(IH1(JL),IH2(JL),IH3(JL)) = XFQSEDH * ZEHS(JL) *  &
                                             ZRHS(JL)**XEXQSEDH * &
                                             ZRHODREFH(JL)**(XEXQSEDH - XCEXVT)
          END IF
        END DO
      END IF
      DO JK = IKB, IKE
        PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1)  - ZWSED(:,:,JK))
        PQHS(:,:,JK) = PQHS(:,:,JK) + ZW(:,:,JK) * (ZWSEDQ(:,:,JK+1) - ZWSEDQ(:,:,JK))
      END DO
      PINPRH(:,:) = PINPRH(:,:) + ZWSED(:,:,IKB)  / XRHOLW / KSPLITR
      IF (JN == KSPLITR) THEN
        PRHS(:,:,:) = PRHS(:,:,:) / PTSTEP
        PQHS(:,:,:) = PQHS(:,:,:) / PTSTEP
      END IF
    END IF
  END DO
!
  IF (OSEDIC) THEN
    IF (ILENALLOCC .GT. 0) DEALLOCATE (ZRCS, ZRHODREFC,  &
        ILISTC,ZWLBDC,ZCONC,ZRCT, ZZT,ZPRES,ZRAY1D,ZFSEDC1D, ZWLBDA,ZCC)         
  END IF
  IF (ILENALLOCR .GT. 0 ) DEALLOCATE(ZRHODREFR,ZRRS,ILISTR)
  IF (ILENALLOCI .GT. 0 ) DEALLOCATE(ZRHODREFI,ZRIS,ILISTI)
  IF (ILENALLOCS .GT. 0 ) DEALLOCATE(ZRHODREFS,ZRSS,ILISTS)
  IF (ILENALLOCG .GT. 0 ) DEALLOCATE(ZRHODREFG,ZRGS,ILISTG)
  IF (KRR == 7 .AND. (ILENALLOCH .GT. 0 )) DEALLOCATE(ZRHODREFH,ZRHS,ILISTH)
!
  IF (ILENALLOCR .GT. 0 ) DEALLOCATE(ZERS,ZQRS,ZLBDAR)
  IF (ILENALLOCI .GT. 0 ) DEALLOCATE(ZEIS,ZQIS,ZCIS,ZCIT,ZLBDAI) 
  IF (ILENALLOCS .GT. 0 ) DEALLOCATE(ZESS,ZQSS,ZLBDAS)
  IF (ILENALLOCG .GT. 0 ) DEALLOCATE(ZEGS,ZQGS,ZLBDAG)
  IF (KRR == 7 .AND. (ILENALLOCH .GT. 0 )) DEALLOCATE(ZEHS,ZQHS,ZLBDAH)
!
!
!*       2.3     budget storage
!
  if ( lbudget_rc .and. osedic ) &
                    call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( osedic ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'SEDI', pqcs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'SEDI', pqrs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'SEDI', pqis(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'SEDI', pqss(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'SEDI', pqgs(:, :, :) * prhodj(:, :, :) )
    if ( krr == 7 ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'SEDI', pqhs(:, :, :) * prhodj(:, :, :) )
  end if
!
  END SUBROUTINE RAIN_ICE_ELEC_SEDIMENTATION_SPLIT
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE RAIN_ICE_SEDIMENTATION_STAT
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!

REAL :: ZP1,ZP2,ZQP,ZH,ZZWLBDA,ZZWLBDC,ZZCC
INTEGER :: JI,JJ,JK
!  
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D !  droplet condensation
!
!------------------------------------------------------------------------------- 
  if ( lbudget_rc .and. osedic ) &
                    call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( osedic ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'SEDI', pqcs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'SEDI', pqrs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'SEDI', pqis(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'SEDI', pqss(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'SEDI', pqgs(:, :, :) * prhodj(:, :, :) )
    if ( krr == 7 ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'SEDI', pqhs(:, :, :) * prhodj(:, :, :) )
  end if
!
!*       1. Parameters for cloud sedimentation
!
  IF (OSEDIC) THEN
    ZRAY(:,:,:)    = 0.
    ZLBC(:,:,:)    = XLBC(1)
    ZFSEDC(:,:,:)  = XFSEDC(1)
    ZCONC3D(:,:,:) = XCONC_LAND
    ZCONC_TMP(:,:) = XCONC_LAND
    IF (PRESENT(PSEA)) THEN
      ZCONC_TMP(:,:) = PSEA(:,:) * XCONC_SEA + (1. - PSEA(:,:)) * XCONC_LAND
      DO JK = IKB, IKE
        ZLBC(:,:,JK)    =  PSEA(:,:) * XLBC(2)   + (1. - PSEA(:,:)) * XLBC(1)
        ZFSEDC(:,:,JK)  = (PSEA(:,:) * XFSEDC(2) + (1. - PSEA(:,:)) * XFSEDC(1))
        ZFSEDC(:,:,JK)  = MAX(MIN(XFSEDC(1),XFSEDC(2)),ZFSEDC(:,:,JK))
        ZCONC3D(:,:,JK) = (1. - PTOWN(:,:)) * ZCONC_TMP(:,:) + PTOWN(:,:) * XCONC_URBAN
        ZRAY(:,:,JK)    = 0.5 * ((1. - PSEA(:,:)) * MOMG(XALPHAC, XNUC, 1.0) + &
                                       PSEA(:,:)  * MOMG(XALPHAC2, XNUC2, 1.0) )
      END DO
    ELSE
      ZCONC3D(:,:,:) = XCONC_LAND                                                     
      ZRAY(:,:,:)    = 0.5 * MOMG(XALPHAC, XNUC, 1.0)
    END IF
    ZRAY(:,:,:)      = MAX(1.,ZRAY(:,:,:))
    ZLBC(:,:,:)      = MAX(MIN(XLBC(1),XLBC(2)),ZLBC(:,:,:))
  ENDIF
!
!
!*       2.    compute the fluxes
!
  ZRTMIN(:) = XRTMIN(:) / PTSTEP
!
  IF (OSEDIC) THEN
    ZPRCS(:,:,:) = 0.0
    ZPRCS(:,:,:) = PRCS(:,:,:) - PRCT(:,:,:) / PTSTEP
    PRCS(:,:,:)  = PRCT(:,:,:) / PTSTEP
  END IF
  ZPRRS(:,:,:) = 0.0
  ZPRSS(:,:,:) = 0.0
  ZPRGS(:,:,:) = 0.0
  IF (KRR == 7) ZPRHS(:,:,:) = 0.0
!
  ZPRRS(:,:,:) = PRRS(:,:,:) - PRRT(:,:,:) / PTSTEP
  ZPRSS(:,:,:) = PRSS(:,:,:) - PRST(:,:,:) / PTSTEP
  ZPRGS(:,:,:) = PRGS(:,:,:) - PRGT(:,:,:) / PTSTEP
  IF (KRR == 7) ZPRHS(:,:,:) = PRHS(:,:,:) - PRHT(:,:,:) / PTSTEP
  PRRS(:,:,:)  = PRRT(:,:,:) / PTSTEP
  PRSS(:,:,:)  = PRST(:,:,:) / PTSTEP
  PRGS(:,:,:)  = PRGT(:,:,:) / PTSTEP
  IF (KRR == 7) PRHS(:,:,:)  = PRHT(:,:,:) / PTSTEP
!
  IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:)
  IF (KRR == 7) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:)
  DO JK = IKB, IKE
    ZW(:,:,JK) = ZTSPLITR / (PRHODREF(:,:,JK) * (PZZ(:,:,JK+1) - PZZ(:,:,JK)))
  END DO
!
!
!*       2.1   for cloud
!
  IF (OSEDIC) THEN
    PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP       
    ZWSED(:,:,:) = 0.
    ZWSEDW1(:,:,:) = 0.
    ZWSEDW2(:,:,:) = 0.

! calculation of P1, P2 and sedimentation flux         
    DO JK = IKE , IKB, -1
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          ! estimation of q' taking into account incomming ZWSED
          ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
          ! calculation of w
          ! mars 2009 : ajout d'un test
          ! IF ( PRCS(JI,JJ,JK) > ZRTMIN(2) ) THEN
          IF(PRCS(JI,JJ,JK) > ZRTMIN(2) .AND. PRCT(JI,JJ,JK) > ZRTMIN(2)) THEN
            ZZWLBDA = 6.6E-8 * (101325. / PPABST(JI,JJ,JK)) * (PTHT(JI,JJ,JK) / 293.15)
            ZZWLBDC = (ZLBC(JI,JJ,JK) * ZCONC3D(JI,JJ,JK) / &
                      (PRHODREF(JI,JJ,JK) * PRCT(JI,JJ,JK)))**XLBEXC
            ZZCC = XCC * (1. + 1.26 * ZZWLBDA * ZZWLBDC / ZRAY(JI,JJ,JK)) ! ZCC: Fall speed
            ZWSEDW1(JI,JJ,JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
                                ZZWLBDC**(-XDC) * ZZCC * ZFSEDC(JI,JJ,JK) 
          ENDIF
          IF (ZQP > ZRTMIN(2)) THEN
            ZZWLBDA = 6.6E-8 * (101325. / PPABST(JI,JJ,JK)) * (PTHT(JI,JJ,JK) / 293.15)
            ZZWLBDC = (ZLBC(JI,JJ,JK) * ZCONC3D(JI,JJ,JK) / &
                      (PRHODREF(JI,JJ,JK) * ZQP))**XLBEXC
            ZZCC = XCC * (1. + 1.26 * ZZWLBDA * ZZWLBDC / ZRAY(JI,JJ,JK)) ! ZCC: Fall speed
            ZWSEDW2(JI,JJ,JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
                                ZZWLBDC**(-XDC) * ZZCC * ZFSEDC(JI,JJ,JK)
          ENDIF
        ENDDO
        DO JI = IIB, IIE
          ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
          ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH)
          ! mars 2009 : correction : ZWSEDW1 =>  ZWSEDW2
          !IF (ZWSEDW1(JI,JJ,JK) /= 0.) THEN
          IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
            ZP2 = MAX(0., 1. -  ZH / (PTSTEP * ZWSEDW2(JI,JJ,JK)) )
          ELSE
            ZP2 = 0.
          ENDIF
          ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) *                          &
                           (PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)) * PRCS(JI,JJ,JK) / &
                            PTSTEP + ZP2 * ZWSED (JI,JJ,JK+1)    
        ENDDO
      ENDDO
    ENDDO  
!
    DO JK = IKB , IKE
      PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
    END DO

    PINPRC(:,:) = ZWSED(:,:,IKB) / XRHOLW         ! in m/s
    PRCS(:,:,:) = PRCS(:,:,:) / PTSTEP
  ENDIF
!
!
!*       2.2   for rain
!
  PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP       
  ZWSED(:,:,:) = 0.
  ZWSEDW1(:,:,:) = 0.
  ZWSEDW2(:,:,:) = 0.
!
! calculation of ZP1, ZP2 and sedimentation flux         
  DO JK = IKE , IKB, -1
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ! estimation of q' taking into account incomming ZWSED
        ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
        ! calculation of w
        IF (PRRS(JI,JJ,JK) > ZRTMIN(3)) THEN
          ZWSEDW1 (JI,JJ,JK) = XFSEDR * PRRS(JI,JJ,JK)**(XEXSEDR-1) * &
                               PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
        ENDIF
        IF (ZQP > ZRTMIN(3)) THEN
          ZWSEDW2(JI,JJ,JK) = XFSEDR * (ZQP)**(XEXSEDR-1) * &
                              PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
        ENDIF
      ENDDO
      DO JI = IIB, IIE
        ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
          ZP2 = MAX(0., 1 -  ZH / (PTSTEP * ZWSEDW2(JI,JJ,JK)) )
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) * &
                          ZH * PRRS(JI,JJ,JK) /      & 
                          PTSTEP + ZP2 * ZWSED (JI,JJ,JK+1)    
      ENDDO
    ENDDO
  ENDDO  

  DO JK = IKB , IKE
    PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
  ENDDO
  PINPRR(:,:) = ZWSED(:,:,IKB) / XRHOLW                        ! in m/s
  PINPRR3D(:,:,:) = ZWSED(:,:,:) / XRHOLW                        ! in m/s
  PRRS(:,:,:) = PRRS(:,:,:) / PTSTEP
!
!
!*       2.3   for pristine ice
!
  PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP       
  ZWSED(:,:,:) = 0.
  ZWSEDW1(:,:,:) = 0.
  ZWSEDW2(:,:,:) = 0.
! calculation of ZP1, ZP2 and sedimentation flux         
  DO JK = IKE , IKB, -1
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ! estimation of q' taking into account incomming ZWSED
        ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
        ! calculation of w
        IF (PRIS(JI,JJ,JK) > MAX(ZRTMIN(4),1.0E-7)) THEN
          ZWSEDW1(JI,JJ,JK) = XFSEDI *                           &
                              PRHODREF(JI,JJ,JK)**(XCEXVT) *     & ! McF&H
                              MAX(0.05E6,-0.15319E6-0.021454E6*  &
                              ALOG(PRHODREF(JI,JJ,JK)*PRIS(JI,JJ,JK)))**XEXCSEDI
        ENDIF
        IF (ZQP > MAX(ZRTMIN(4),1.0E-7)) THEN
          ZWSEDW2(JI,JJ,JK)= XFSEDI *                           &
                             PRHODREF(JI,JJ,JK)**(XCEXVT) *     & ! McF&H
                             MAX( 0.05E6,-0.15319E6-0.021454E6* &
                             ALOG(PRHODREF(JI,JJ,JK)*ZQP) )**XEXCSEDI
        ENDIF
      ENDDO
      DO JI = IIB, IIE
        ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH / (PTSTEP * ZWSEDW2(JI,JJ,JK)))
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) *                          &
                         (PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)) * PRIS(JI,JJ,JK) / &
                          PTSTEP + ZP2 * ZWSED(JI,JJ,JK+1)    
      ENDDO
    ENDDO
  ENDDO  
!
  DO JK = IKB , IKE
    PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
  ENDDO
!
  PRIS(:,:,:) = PRIS(:,:,:) / PTSTEP
!
!
!*       2.4   for aggregates/snow
!
  PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP       
  ZWSED(:,:,:) = 0.
  ZWSEDW1(:,:,:) = 0.
  ZWSEDW2(:,:,:) = 0.

! calculation of ZP1, ZP2 and sedimentation flux         
  DO JK = IKE , IKB, -1
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ! estimation of q' taking into account incomming ZWSED 
         ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
        ! calculation of w
        IF (PRSS(JI,JJ,JK) > ZRTMIN(5)) THEN
          ZWSEDW1(JI,JJ,JK) = XFSEDS * (PRSS(JI,JJ,JK))**(XEXSEDS-1) * &
                              PRHODREF(JI,JJ,JK)**(XEXSEDS-XCEXVT-1)
        ENDIF
        IF (ZQP > ZRTMIN(5)) THEN
          ZWSEDW2(JI,JJ,JK) = XFSEDS * (ZQP)**(XEXSEDS-1) * &
                              PRHODREF(JI,JJ,JK)**(XEXSEDS-XCEXVT-1)
        ENDIF
      ENDDO
      DO JI = IIB, IIE
        ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH / (PTSTEP * ZWSEDW2(JI,JJ,JK)) )
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) * &
                          ZH * PRSS(JI,JJ,JK) /      &
                          PTSTEP + ZP2 * ZWSED(JI,JJ,JK+1)    
      ENDDO
    ENDDO
  ENDDO  
!
  DO JK = IKB , IKE
    PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
  ENDDO
!
  PINPRS(:,:) = ZWSED(:,:,IKB) / XRHOLW      ! in m/s
  PRSS(:,:,:) = PRSS(:,:,:) / PTSTEP
!
!
!
!*       2.5   for graupeln
!
  PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP       
  ZWSED(:,:,:) = 0.
  ZWSEDW1(:,:,:) = 0.
  ZWSEDW2(:,:,:) = 0.

! calculation of ZP1, ZP2 and sedimentation flux         
  DO JK = IKE , IKB, -1
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ! estimation of q' taking into account incomming ZWSED 
        ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
        ! calculation of w
        IF (PRGS(JI,JJ,JK) > ZRTMIN(6)) THEN
          ZWSEDW1(JI,JJ,JK) = XFSEDG * (PRGS(JI,JJ,JK))**(XEXSEDG-1) * &
                              PRHODREF(JI,JJ,JK)**(XEXSEDG-XCEXVT-1)
        ENDIF
        IF (ZQP > ZRTMIN(6)) THEN
            ZWSEDW2(JI,JJ,JK) = XFSEDG * (ZQP)**(XEXSEDG-1) * &
                                PRHODREF(JI,JJ,JK)**(XEXSEDG-XCEXVT-1)
        ENDIF
      ENDDO
      DO JI = IIB, IIE
        ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH / (PTSTEP * ZWSEDW2(JI,JJ,JK)) )
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) * &
                          ZH * PRGS(JI,JJ,JK) /      &
                          PTSTEP + ZP2 * ZWSED(JI,JJ,JK+1)    
      ENDDO
    ENDDO
  ENDDO  
!
  DO JK = IKB , IKE
    PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
  ENDDO

  PINPRG(:,:) = ZWSED(:,:,IKB) / XRHOLW       ! in m/s
  PRGS(:,:,:) = PRGS(:,:,:) / PTSTEP
!
!
!*       2.6   for hail
!
  IF (KRR == 7) THEN
    PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP       
    ZWSED(:,:,:) = 0.
    ZWSEDW1(:,:,:) = 0.
    ZWSEDW2(:,:,:) = 0.
! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = IKE , IKB, -1
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          ! estimation of q' taking into account incomming ZWSED 
          ZQP = ZWSED(JI,JJ,JK+1) * ZW(JI,JJ,JK)
          ! calculation of w
          IF ((PRHS(JI,JJ,JK)+ZQP) > ZRTMIN(7) ) THEN
            ZWSEDW1 (JI,JJ,JK) = XFSEDH * (PRHS(JI,JJ,JK))**(XEXSEDH-1) * &
                                 PRHODREF(JI,JJ,JK)**(XEXSEDH-XCEXVT-1)
        ENDIF
        IF (ZQP > ZRTMIN(7)) THEN
          ZWSEDW2(JI,JJ,JK) = XFSEDH * ZQP**(XEXSEDH-1) * &
                              PRHODREF(JI,JJ,JK)**(XEXSEDH-XCEXVT-1)
          ENDIF
        ENDDO
        DO JI = IIB, IIE
          ZH = PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
          ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH)
          IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
            ZP2 = MAX(0.,1 - ZH / (PTSTEP*ZWSEDW2(JI,JJ,JK)))
          ELSE
            ZP2 = 0.
          ENDIF
          ZWSED(JI,JJ,JK) = ZP1 * PRHODREF(JI,JJ,JK) * &
                            ZH * PRHS(JI,JJ,JK) /      &
                            PTSTEP + ZP2 * ZWSED(JI,JJ,JK+1)
        ENDDO
      ENDDO
    ENDDO
!
    DO JK = IKB , IKE
      PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK) * (ZWSED(:,:,JK+1) - ZWSED(:,:,JK))
    ENDDO
!
    PINPRH(:,:) = ZWSED(:,:,IKB) / XRHOLW       ! in m/s
    PRHS(:,:,:) = PRHS(:,:,:) / PTSTEP
  ENDIF
!
!
!*       2.3     budget storage
!
  if ( lbudget_rc .and. osedic ) &
                    call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( osedic ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'SEDI', pqcs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'SEDI', pqrs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'SEDI', pqis(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'SEDI', pqss(:, :, :) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'SEDI', pqgs(:, :, :) * prhodj(:, :, :) )
    if ( krr == 7 ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'SEDI', pqhs(:, :, :) * prhodj(:, :, :) )
  end if
!
  END SUBROUTINE RAIN_ICE_SEDIMENTATION_STAT
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_ELEC_NUCLEATION
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
INTEGER , DIMENSION(SIZE(GNEGT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
!  compute the temperature and the pressure
!
ZT(:,:,:) = PTHT(:,:,:) * (PPABST(:,:,:) / XP00) ** (XRD / XCPD)
!
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:) = .FALSE.
GNEGT(IIB:IIE,IJB:IJE,IKB:IKE) = ZT(IIB:IIE,IJB:IJE,IKB:IKE) < XTT
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
!
IF( INEGT >= 1 ) THEN
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HENU', pris(:, :, :) * prhodj(:, :, :) )

  ALLOCATE(ZRVT(INEGT))
  ALLOCATE(ZCIT(INEGT)) 
  ALLOCATE(ZZT(INEGT))  
  ALLOCATE(ZPRES(INEGT))
  DO JL = 1, INEGT
    ZRVT(JL)  = PRVT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL)  = PCIT(I1(JL),I2(JL),I3(JL))
    ZZT(JL)   = ZT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
  ENDDO
  ALLOCATE(ZZW(INEGT))
  ALLOCATE(ZUSW(INEGT))
  ALLOCATE(ZSSI(INEGT))
  ZZW(:) = EXP(XALPI - XBETAI / ZZT(:) - XGAMI * ALOG(ZZT(:)))         ! es_i
  ZZW(:) = MIN(ZPRES(:) / 2., ZZW(:))             ! safety limitation
  ZSSI(:) = ZRVT(:) * (ZPRES(:) - ZZW(:)) / ((XMV / XMD) * ZZW(:)) - 1.0
                                                  ! Supersaturation over ice
  ZUSW(:) = EXP(XALPW - XBETAW / ZZT(:) - XGAMW * ALOG(ZZT(:)))        ! es_w
  ZUSW(:) = MIN(ZPRES(:) / 2., ZUSW(:))           ! safety limitation
  ZUSW(:) = (ZUSW(:) / ZZW(:)) * ((ZPRES(:) - ZZW(:)) / (ZPRES(:) - ZUSW(:))) - 1.0
                           ! Supersaturation of saturated water vapor over ice
!
!*       3.1     compute the heterogeneous nucleation source: RVHENI
!
!*       3.1.1   compute the cloud ice concentration
!
  ZZW(:) = 0.0
  ZSSI(:) = MIN( ZSSI(:), ZUSW(:) ) ! limitation of SSi according to SSw=0
!
  WHERE ((ZZT(:) < XTT-5.0) .AND. (ZSSI(:) > 0.0))
    ZZW(:) = XNU20 * EXP(XALPHA2 * ZSSI(:) - XBETA2)
  END WHERE
  WHERE ((ZZT(:) <= XTT-2.0) .AND. (ZZT(:) >= XTT-5.0) .AND. (ZSSI(:) > 0.0))
    ZZW(:) = MAX(XNU20 * EXP(-XBETA2), XNU10 * EXP(-XBETA1 * (ZZT(:) - XTT)) * &
                               (ZSSI(:) / ZUSW(:))**XALPHA1 )
  END WHERE
  ZZW(:) = ZZW(:) - ZCIT(:)
!
  IF( MAXVAL(ZZW(:)) > 0.0 ) THEN
!
!*       3.1.2   update the r_i and r_v mixing ratios
!
    ZZW(:) = MIN(ZZW(:), 50.E3) ! limitation provisoire a 50 l^-1
    ZW(:,:,:) = UNPACK(ZZW(:), MASK=GNEGT(:,:,:), FIELD=0.0)
    ZW(:,:,:) = MAX(ZW(:,:,:), 0.0) * XMNU0 / (PRHODREF(:,:,:) * PTSTEP)
    PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
    PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
    IF (KRR == 7) THEN
      PTHS(:,:,:) = PTHS(:,:,:) +                                                    &
                    ZW(:,:,:) * (XLSTT + (XCPV - XCI) * (ZT(:,:,:) - XTT)) /         &
                  ((XCPD + XCPV * PRVT(:,:,:) + XCL * (PRCT(:,:,:) + PRRT(:,:,:)) +  &
                    XCI * (PRIT(:,:,:) + PRST(:,:,:) + PRGT(:,:,:) + PRHT(:,:,:))) * &
                    PEXNREF(:,:,:))
    ELSE IF(KRR == 6) THEN
      PTHS(:,:,:) = PTHS(:,:,:) + &
                    ZW(:,:,:) * (XLSTT + (XCPV - XCI) * (ZT(:,:,:) - XTT))  /       &
                  ((XCPD + XCPV * PRVT(:,:,:) + XCL * (PRCT(:,:,:) + PRRT(:,:,:)) + &
                    XCI * (PRIT(:,:,:) + PRST(:,:,:) + PRGT(:,:,:))) * PEXNREF(:,:,:))
    END IF
! f(L_s*(RVHENI))
    ZZW(:) = MAX( ZZW(:)+ZCIT(:),ZCIT(:) )
    PCIT(:,:,:) = MAX( UNPACK( ZZW(:),MASK=GNEGT(:,:,:),FIELD=0.0 ) , &
                       PCIT(:,:,:) )
  END IF
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZRVT)

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HENU', pris(:, :, :) * prhodj(:, :, :) )

END IF

  END SUBROUTINE RAIN_ICE_ELEC_NUCLEATION
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE RAIN_ICE_ELEC_SLOW
!
!*      0. DECLARATIONS
!          ------------
USE MODD_CST, ONLY : XMNH_HUGE_12_LOG
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!
!*       3.5.1   compute the homogeneous nucleation source: RCHONI & QCHONI
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'HON', &
                            Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'HON', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
  ZWQ1(:,1:7) = 0.0
!
  WHERE( ABS(ZECT(:)) <= XECMIN)
    ZECT(:) = 0.
  ENDWHERE
!
  WHERE( (ZZT(:)<XTT-35.0) .AND. (ZRCT(:)>XRTMIN(2)) .AND. (ZRCS(:)>0.) )
    ZZW(:) = MIN( ZRCS(:),XHON*ZRHODREF(:)*ZRCT(:)       &
                                 *EXP( MIN(XMNH_HUGE_12_LOG,XALPHA3*(ZZT(:)-XTT)-XBETA3) ) )
    ZRIS(:) = ZRIS(:) + ZZW(:)
    ZRCS(:) = ZRCS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCHONI))
    ZWQ1(:,1) = XQHON * ZECT(:) * ZZW(:)  ! QCHONI
  ENDWHERE
!
  WHERE (ZZT(:) < (XTT - 35.)     .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRCS(:) > ZRSMIN_ELEC(2) .AND. &
         ABS(ZQCS(:)) > XQTMIN(2) .AND. ABS(ZECT(:)) > XECMIN)                    
    ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
    ZQIS(:) = ZQIS(:) + ZWQ1(:,1)
    ZQCS(:) = ZQCS(:) - ZWQ1(:,1) 
  END WHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HON', Unpack(  zzw(:) * ( zlsfact(:) - zlvfact(:) ) &
                                                           * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'HON', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'HON', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'HON', &
                           Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'HON', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
  end if
!
!*       3.5.2   compute the spontaneous freezing source: RRHONG & QRHONG
!
  ZZW(:) = 0.0
!
  WHERE( (ZZT(:)<XTT-35.0) .AND. (ZRRT(:)>XRTMIN(3)) .AND. (ZRRS(:)>0.) )
    ZZW(:) = MIN( ZRRS(:),ZRRT(:)/PTSTEP )
    ZRGS(:) = ZRGS(:) + ZZW(:)
    ZRRS(:) = ZRRS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRHONG))
  ENDWHERE
!
  WHERE (ZZT(:) < (XTT - 35.) .AND.                                    &
         ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
         ZRGS(:) > ZRSMIN_ELEC(6) .AND. ABS(ZQRT(:)) > XQTMIN(3))
    ZWQ1(:,2) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZQRT(:)/PTSTEP) ),ZQRS(:) ) ! QRHONG
    ZQGS(:) = ZQGS(:) + ZWQ1(:,2)
    ZQRS(:) = ZQRS(:) - ZWQ1(:,2)
  ENDWHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'SFR', Unpack(  zzw(:) * ( zlsfact(:) - zlvfact(:) ) &
                                                           * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'SFR', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'SFR', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'SFR', &
                           Unpack( -zwq1(:, 2) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'SFR', &
                           Unpack(  zwq1(:, 2) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
!
!*       3.5.3  compute the deposition, aggregation and autoconversion sources
!
  ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
  ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
!
!*       3.5.3.1  compute the thermodynamical function A_i(T,P)
!*                and the c^prime_j (in the ventilation factor)
!
  ZAI(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
  ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                 + ( XRV*ZZT(:) ) / (ZDV(:)*ZAI(:))
  ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!*      3.5.3.2  compute the riming-conversion of r_c for r_i production: RCAUTI
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'DEPS', &
                            Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'DEPS', &
                            Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'DEPS', &
                            Unpack( zqss(:) * zrhodj(:),  mask = gmicro(:, :, :), field =  0. ) )
  end if

  ZZW(:) = 0.0
!
  WHERE ((ZRST(:) > XRTMIN(5)) .AND. (ZRSS(:) > 0.0))
    ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *                               &
             ( X0DEPS*ZLBDAS(:)**XEX0DEPS + X1DEPS*ZCJ(:)*ZLBDAS(:)**XEX1DEPS )
    ZZW(:) =         MIN( ZRVS(:),ZZW(:)      ) * (0.5 + SIGN(0.5,ZZW(:))) &
                   - MIN( ZRSS(:),ABS(ZZW(:)) ) * (0.5 - SIGN(0.5,ZZW(:)))
    ZRSS(:) = ZRSS(:) + ZZW(:)
    ZRVS(:) = ZRVS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
    ZWQ1(:,5) = XCOEF_RQ_S * ZQST(:) * (-ZZW(:)) / ZRST(:) ! sublimation
  END WHERE
!
  WHERE (ZRST(:) > XRTMIN_ELEC(5) .AND. ZRSS(:) > ZRSMIN_ELEC(5) .AND. &
         ZRVS(:) > ZRSMIN_ELEC(1) .AND. ABS(ZQST(:)) > XQTMIN(5) .AND. &
         ZZW(:) < 0. .AND. (-ZZW(:) <= ZRSS(:)))
    ZWQ1(:,5) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,5)) ),ZQSS(:) )
    ZQSS(:) = ZQSS(:) - ZWQ1(:,5)
    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ1(:,5)/XECHARGE )
    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ1(:,5)/XECHARGE )
  ENDWHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPS', Unpack(  zzw(:) * zlsfact(:) &
                                                           * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPS', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'DEPS', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'DEPS', &
                           Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'DEPS', &
                           Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'DEPS', &
                           Unpack( zqss(:) * zrhodj(:),  mask = gmicro(:, :, :), field =  0. ) )
  end if
!
!*       3.5.3.4  compute the aggregation on r_s: RIAGGS & QIAGGS
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'AGGS', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'AGGS', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
  end if

  ZZW(:) = 0.0
  ZLATHAMIAGGS(:) = 1.0
  IF (LIAGGS_LATHAM) THEN
    ZLATHAMIAGGS(:) = 1.0 + 0.4E-10 * MIN( 2.25E10,                &
                      ZEFIELDU(:)**2+ZEFIELDV(:)**2+ZEFIELDW(:)**2 )
  ENDIF
!
  WHERE (ZRIT(:) > XRTMIN(4) .AND. ZRST(:) > XRTMIN(5) .AND. ZRIS(:) > 0.0)
    ZZW(:) = MIN( ZRIS(:),XFIAGGS * EXP( XCOLEXIS*(ZZT(:)-XTT) ) &
                    * ZLATHAMIAGGS(:)                            &
                    * ZRIT(:)                                    &
                    * ZLBDAS(:)**XEXIAGGS                        &
                    * ZRHOCOR(:) / ZCOR00                        )
    ZRSS(:)  = ZRSS(:) + ZZW(:)
    ZRIS(:)  = ZRIS(:) - ZZW(:)
    ZWQ1(:,3) = XCOEF_RQ_I * ZZW(:) * ZQIT(:) / ZRIT(:) ! QIAGGS_coal
  END WHERE
!
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
         ZRSS(:) > ZRSMIN_ELEC(5) .AND. ABS(ZQIT(:)) > XQTMIN(4))
    ZWQ1(:,3) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,3)) ),ZQIS(:) )
    ZQSS(:) = ZQSS(:) + ZWQ1(:,3)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,3)
  END WHERE

  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AGGS', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AGGS', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'AGGS', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'AGGS', &
                           Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
  end if

  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'NIIS', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'NIIS', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field =  0. ) )
  end if

  CALL ELEC_IAGGS_B()                                   ! QIAGGS_boun

  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'NIIS', &
                           Unpack( zqis(:), mask = gmicro(:, :, :), field =  0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'NIIS', &
                           Unpack( zqss(:), mask = gmicro(:, :, :), field =  0. ) )
  end if

! Save the NI charging rate for temporal series
  XNI_IAGGS(:,:,:) = UNPACK(ZWQ1(:,7), MASK=GMICRO, FIELD=0.0)
  XNI_IAGGS(:,:,:) = XNI_IAGGS(:,:,:) * PRHODREF(:,:,:)  ! C/m3/s
!
!*       3.5.3.5  compute the autoconversion of r_i for r_s production: 
!                 RIAUTS & QIAUTS
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'AUTS', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'AUTS', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ALLOCATE( ZCRIAUTI(IMICRO ))
  ZCRIAUTI(:) = MIN(XCRIAUTI,10**(0.06*(ZZT(:)-XTT)-3.5))
  ZZW(:) = 0.0
!
  WHERE ((ZRIT(:) > XRTMIN(4)) .AND. (ZRIS(:) > 0.0))
    ZZW(:) = MIN( ZRIS(:),XTIMAUTI * EXP( XTEXAUTI*(ZZT(:)-XTT) ) &
                            * MAX( ZRIT(:)-ZCRIAUTI(:),0.0 ) )
    ZRSS(:) = ZRSS(:) + ZZW(:)
    ZRIS(:) = ZRIS(:) - ZZW(:)
    ZWQ1(:,4) = XCOEF_RQ_I * ZQIT(:) * ZZW(:) / ZRIT(:)  ! QIAUTS
  END WHERE
!
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
         ZRSS(:) > ZRSMIN_ELEC(5) .AND. ABS(ZQIT(:)) > XQTMIN(4))
    ZWQ1(:,4) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,4)) ),ZQIS(:) )
    ZQSS(:) = ZQSS(:) + ZWQ1(:,4)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,4)
  END WHERE
!
  DEALLOCATE(ZCRIAUTI)

  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AUTS', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AUTS', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'AUTS', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'AUTS', &
                           Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
!
!*       3.5.3.6  compute the deposition on r_g: RVDEPG & QVDEPG
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'DEPG', &
                            Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'DEPG', &
                            Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'DEPG', &
                            Unpack( zqgs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
!
  WHERE ((ZRGT(:) > XRTMIN(6)) .AND. (ZRGS(:) > 0.0))
    ZZW(:) = (ZSSI(:) / (ZRHODREF(:) * ZAI(:))) *                              &
             (X0DEPG * ZLBDAG(:)**XEX0DEPG + X1DEPG * ZCJ(:) * ZLBDAG(:)**XEX1DEPG)
    ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                   - MIN( ZRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
    ZRGS(:) = ZRGS(:) + ZZW(:)
    ZRVS(:) = ZRVS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
    ZWQ1(:,6) = XCOEF_RQ_G * ZQGT(:) * (-ZZW(:)) / ZRGT(:)      ! sublimation
  END WHERE
!
  WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND. &
         ZRVS(:) > ZRSMIN_ELEC(1) .AND. ABS(ZQGT(:)) > XQTMIN(6) .AND. &
         ZZW(:) < 0. .AND. (-ZZW(:)) <= ZRGS(:))
    ZWQ1(:,6) = SIGN( MIN( ABS(ZQGS(:)),ABS(ZWQ1(:,6)) ),ZQGS(:) )
    ZQGS(:) = ZQGS(:) - ZWQ1(:,6) 
    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ1(:,6)/XECHARGE )
    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ1(:,6)/XECHARGE )
  END WHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPG', Unpack(  zzw(:) * zlsfact(:) &
                                                           * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPG', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DEPG', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'DEPG', &
                           Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'DEPG', &
                           Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'DEPG', &
                           Unpack( zqgs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
  end if

  END SUBROUTINE RAIN_ICE_ELEC_SLOW
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_ELEC_WARM
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
REAL :: ZCRIAUTC             ! Critical cloud mixing ratio
!
!-------------------------------------------------------------------------------
!
!*       4.1    compute the autoconversion of r_c for r_r production: 
!               RCAUTR & QCAUTR
!
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'AUTO', &
                              Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'AUTO', &
                              Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if

    ZZW(:) = 0.0
    ZWQ1(:,1:3) = 0.0
!
    IF ( HSUBG_AUCV == 'CLFR' ) THEN
      WHERE ((ZRCT(:) > 0.0) .AND. (ZRCS(:) > 0.0) .AND. (ZCF(:) > 0.0))
        ZZW(:) = XTIMAUTC * MAX( ZRCT(:)/(ZCF(:)) -XCRIAUTC/ZRHODREF(:),0.0)
        ZZW(:) = MIN( ZRCS(:),(ZCF(:))*ZZW(:))
        ZRCS(:) = ZRCS(:) - ZZW(:)
        ZRRS(:) = ZRRS(:) + ZZW(:)
        ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW(:) / ZRCT(:)   ! QCAUTR
      END WHERE
    ELSE IF (HSUBG_AUCV == 'SIGM') THEN
      DO JL = 1, IMICRO
        IF (ZRCS(JL) > 0.0) THEN
          ZCRIAUTC = XCRIAUTC / ZRHODREF(JL)
          IF (ZRCT(JL) > (ZCRIAUTC + ZSIGMA_RC(JL))) THEN
            ZZW(JL) = MIN( ZRCS(JL) , XTIMAUTC* ( ZRCT(JL)-ZCRIAUTC ) )
          ELSEIF (ZRCT(JL) >  (ZCRIAUTC - ZSIGMA_RC(JL)) .AND. &
                  ZRCT(JL) <= (ZCRIAUTC + ZSIGMA_RC(JL))) THEN
            ZZW(JL) = MIN( ZRCS(JL) , XTIMAUTC*( ZRCT(JL)+ZSIGMA_RC(JL)-ZCRIAUTC )**2 &
                                                /( 4. * ZSIGMA_RC(JL) )                 )
          ENDIF
          ZRCS(JL) = ZRCS(JL) - ZZW(JL)
          ZRRS(JL) = ZRRS(JL) + ZZW(JL)
          IF (ZRCT(JL) > 0.) THEN
            ZWQ1(JL,1) = XCOEF_RQ_C * ZQCT(JL) * ZZW(JL) / ZRCT(JL) 
          END IF
        ENDIF
      END DO
    ELSE
      WHERE ((ZRCT(:) > XRTMIN(2)) .AND. (ZRCS(:) > 0.0))
        ZZW(:) = MIN( ZRCS(:),XTIMAUTC*MAX( ZRCT(:)-XCRIAUTC/ZRHODREF(:),0.0 ) )
        ZRCS(:) = ZRCS(:) - ZZW(:)
        ZRRS(:) = ZRRS(:) + ZZW(:)
        ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW(:) / ZRCT(:)     ! QCAUTR
      END WHERE
    END IF
!
    WHERE (ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRCS(:) > ZRSMIN_ELEC(2) .AND. &
           ZRRS(:) > ZRSMIN_ELEC(3) .AND. ABS(ZQCT(:)) > XQTMIN(2))
      ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
      ZQCS(:) = ZQCS(:) - ZWQ1(:,1)
      ZQRS(:) = ZQRS(:) + ZWQ1(:,1)
    END WHERE

    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'AUTO', &
                                             Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'AUTO', &
                                             Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'AUTO', &
                             Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'AUTO', &
                             Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
!
!*       4.2    compute the accretion of r_c for r_r production: RCACCR & QCACCR
!
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'ACCR', &
                              Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'ACCR', &
                              Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if

    ZZW(:) = 0.0
    WHERE ((ZRCT(:) > XRTMIN(2)) .AND. (ZRRT(:) > XRTMIN(3)) .AND. (ZRCS(:) > 0.0))
      ZZW(:) = MIN( ZRCS(:),XFCACCR * ZRCT(:)             &
                                    * ZLBDAR(:)**XEXCACCR &
                                    * ZRHOCOR(:)/ZCOR00 )
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRRS(:) = ZRRS(:) + ZZW(:)
      ZWQ1(:,2) = XCOEF_RQ_C * ZQCT(:) * ZZW(:) / ZRCT(:)          ! QCACCR
    END WHERE
!
    WHERE (ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
           ZRCS(:) > ZRSMIN_ELEC(2) .AND. ABS(ZQCT(:)) > XQTMIN(2))
      ZWQ1(:,2) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,2)) ),ZQCS(:) )
      ZQCS(:) = ZQCS(:) - ZWQ1(:,2)
      ZQRS(:) = ZQRS(:) + ZWQ1(:,2)
    ENDWHERE

    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'ACCR', &
                                             Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'ACCR', &
                                             Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'ACCR', &
                             Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'ACCR', &
                             Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
!
!
!*       4.3    compute the evaporation of r_r: RREVAV & QREVAV
!
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'REVA', &
                              Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'REVA', &
                              Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'REVA', &
                              Unpack( zqrs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
    end if

    ZZW(:) = 0.0
    WHERE ((ZRRT(:) > XRTMIN(3)) .AND. (ZRCT(:) <= XRTMIN(2)))
      ZZW(:)  = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
      ZUSW(:) = 1.0 - ZRVT(:) * (ZPRES(:) - ZZW(:)) / ((XMV / XMD) * ZZW(:)) 
                                                    ! Undersaturation over water
      ZZW(:) = (XLVTT + (XCPV - XCL) * (ZZT(:) - XTT) )**2 / &
               (ZKA(:) * XRV * ZZT(:)**2) +                  &
               (XRV * ZZT(:)) / (ZDV(:) * ZZW(:))
      ZZW(:) = MIN( ZRRS(:),( MAX( 0.0,ZUSW(:) )/(ZRHODREF(:)*ZZW(:)) ) *      &
              ( X0EVAR*ZLBDAR(:)**XEX0EVAR+X1EVAR*ZCJ(:)*ZLBDAR(:)**XEX1EVAR ) )
      ZRRS(:) = ZRRS(:) - ZZW(:)
      ZRVS(:) = ZRVS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) - ZZW(:)*ZLVFACT(:)
      ZWQ1(:,3) = XCOEF_RQ_R * ZQRT(:) * ZZW(:) / ZRRT(:)     ! QREVAV
    END WHERE
!
    WHERE (ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
           ZRVS(:) > ZRSMIN_ELEC(1) .AND. ZRCT(:) <= 0.0           .AND. &
           ABS(ZQRT(:)) > XQTMIN(3)) 
      ZWQ1(:,3) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,3)) ),ZQRS(:) )
      ZQRS(:) = ZQRS(:) - ZWQ1(:,3)
      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ1(:,3)/XECHARGE )
      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ1(:,3)/XECHARGE )
    ENDWHERE
!
    PEVAP3D(:,:,:)=UNPACK(ZZW(:),MASK=GMICRO(:,:,:),FIELD=PEVAP3D(:,:,:))

    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'REVA', &
                                Unpack( -zzw(:) * zlvfact(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'REVA', &
                                             Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'REVA', &
                                             Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ), 'REVA', &
                             Unpack( zqpis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend ), 'REVA', &
                             Unpack( zqnis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'REVA', &
                             Unpack( zqrs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
    end if

  END SUBROUTINE RAIN_ICE_ELEC_WARM
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE RAIN_ICE_ELEC_FAST_RS
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!
  ZZW1(:,:)   = 0.0
  ZWQ1(:,1:7) = 0.0
!
  ALLOCATE( GRIM(IMICRO) )
  GRIM(:) = (ZRCT(:) > XRTMIN(2)) .AND. (ZRST(:) > XRTMIN(5)) .AND. &
            (ZRCS(:) > 0.0) .AND. (ZZT(:) < XTT)
  IGRIM = COUNT( GRIM(:) )
!
  IF (IGRIM > 0) THEN
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'RIM', &
                              Unpack( zqcs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'RIM', &
                              Unpack( zqss(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'RIM', &
                              Unpack( zqgs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
    end if
!
!        5.1.0  allocations
!
    ALLOCATE( ZVEC1(IGRIM) )
    ALLOCATE( ZVEC2(IGRIM) )
    ALLOCATE( IVEC1(IGRIM) )
    ALLOCATE( IVEC2(IGRIM) )
!
!*       5.1.1  select the ZLBDAS
!
    ZVEC1(:) = PACK( ZLBDAS(:),MASK=GRIM(:) )
!
!*       5.1.2  find the next lower indice for the ZLBDAS in the geometrical
!*              set of Lbda_s used to tabulate some moments of the incomplete 
!               gamma function
!
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(NGAMINC)-0.00001,           &
                          XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
!
!*       5.1.3  perform the linear interpolation of the normalized
!*              "2+XDS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!*       5.1.4  riming of the small sized aggregates
!
    WHERE (GRIM(:) .AND. ZRCS(:) > 0.0) 
      ZZW1(:,1) = MIN( ZRCS(:),                      &
                       XCRIMSS * ZZW(:) * ZRCT(:) *  & ! RCRIMSS
                       ZLBDAS(:)**XEXCRIMSS * ZRHOCOR(:)/ZCOR00 )
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRSS(:) = ZRSS(:) + ZZW1(:,1)
      ZTHS(:) = ZTHS(:) + ZZW1(:,1) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*(RCRIMSS))
      ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW1(:,1) / ZRCT(:)     ! QCRIMSS
    END WHERE
!
    WHERE (ZZT(:) < XTT .AND.                                            &
           ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRSS(:) > ZRSMIN_ELEC(5) .AND. &
           ABS(ZQCT(:)) > XQTMIN(2) .AND. ZRCS(:) > ZRSMIN_ELEC(2)) 
      ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
      ZQCS(:) = ZQCS(:) - ZWQ1(:,1)
      ZQSS(:) = ZQSS(:) + ZWQ1(:,1)
    ENDWHERE
!
!*       5.1.5  perform the linear interpolation of the normalized
!*              "XBS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!
!*	 5.1.6	perform the linear interpolation of the normalized
!*		"XFS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =  XGAMINC_RIM3( IVEC2(1:IGRIM)+1 ) *  ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM3( IVEC2(1:IGRIM)   ) * (ZVEC2(1:IGRIM) - 1.0)
    ZWQ1(:,3) = UNPACK( VECTOR=ZVEC1(:), MASK=GRIM, FIELD=0.0 )
!
!*       5.1.7  riming-conversion of the large sized aggregates into graupeln:
!*              RSRIMCG & QSRIMCG and RCRIMSG & QCRIMSG
!
    WHERE (GRIM(:) .AND. ZRSS(:) > 0.0 .AND. ZRCS(:) > 0.0 .AND. ZZW(:) < 1.) 
      ZZW1(:,2) = MIN( ZRCS(:),                         &
                       XCRIMSG * ZRCT(:)                & ! RCRIMSG
                               * ZLBDAS(:)**XEXCRIMSG   &
                               * ZRHOCOR(:)/ZCOR00 - ZZW1(:,1) )
      ZZW1(:,3) = MIN( ZRSS(:),                         &
                       XSRIMCG * ZLBDAS(:)**XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(:) )/(PTSTEP*ZRHODREF(:)) )
      ZRCS(:) = ZRCS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2) + ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*(RCRIMSG))
      ZWQ1(:,2) = XCOEF_RQ_C * ZQCT(:) * ZZW1(:,2) / ZRCT(:)      ! QCRIMSG
      ZWQ1(:,3) = XQSRIMCG * ZEST(:) *                       &    ! QSRIMCG
                  ZLBDAS(:)**XEXQSRIMCG * (1. - ZWQ1(:,3)) / &
                  (PTSTEP * ZRHODREF(:))
    END WHERE
!
    WHERE (ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZRCS(:) > ZRSMIN_ELEC(2) .AND. &
           ZZT(:) < XTT             .AND. ABS(ZQCT(:)) > XQTMIN(2))
      ZWQ1(:,2) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,2)) ),ZQCS(:) )
      ZQGS(:) = ZQGS(:) + ZWQ1(:,2)
      ZQCS(:) = ZQCS(:) - ZWQ1(:,2)
    ENDWHERE
!
    WHERE (ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZRCS(:) > ZRSMIN_ELEC(2) .AND. &
           ZZT(:) < XTT             .AND. ABS(ZQCT(:)) > XQTMIN(2) .AND. &
           ABS(ZEST) > XESMIN)
      ZWQ1(:,3) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,3)) ),ZQSS(:) )
      ZQGS(:) = ZQGS(:) + ZWQ1(:,3) 
      ZQSS(:) = ZQSS(:) - ZWQ1(:,3)
    ENDWHERE
!
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)

    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'RIM', Unpack( ( zzw1(:,1) + zzw1(:,2) ) &
                                                  * ( zlsfact(:) - zlvfact(:) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'RIM', &
                                             Unpack( ( -zzw1(:,1) - zzw1(:,2) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'RIM', &
                                             Unpack( (  zzw1(:,1) - zzw1(:,3) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'RIM', &
                                             Unpack( (  zzw1(:,2) + zzw1(:,3) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'RIM', &
                             Unpack( zqcs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'RIM', &
                             Unpack( zqss(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'RIM', &
                             Unpack( zqgs(:) * zrhodj(:),  mask = gmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
  DEALLOCATE(GRIM)
!
!
!*       5.2    rain accretion onto the aggregates
!
  ZZW1(:,2:3) = 0.0
  ZWQ4(:)     = 0.0
!
  ALLOCATE(GACC(IMICRO))
  GACC(:) = ZRRT(:)>XRTMIN(3) .AND. ZRST(:)>XRTMIN(5) .AND. &
            ZRRS(:) > 0.0 .AND. ZZT(:) < XTT
  IGACC = COUNT( GACC(:) )
!
  IF( IGACC>0 ) THEN
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'ACC', &
                                              Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACC', &
                                              Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'ACC', &
                                              Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'ACC', &
                                              Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'ACC', &
                              Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'ACC', &
                              Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'ACC', &
                              Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
!
!        5.2.0  allocations
!
    ALLOCATE(ZVEC1(IGACC))
    ALLOCATE(ZVEC2(IGACC))
    ALLOCATE(ZVEC3(IGACC))
    ALLOCATE(IVEC1(IGACC))
    ALLOCATE(IVEC2(IGACC))
!
    ALLOCATE( ZVECQ4(IGACC) )
    ALLOCATE( ZVECQ5(IGACC) )
    ALLOCATE( ZVECQ6(IGACC) )
!
!
!        5.2.1  select the (ZLBDAS,ZLBDAR) couplet
!
    ZVEC1(:) = PACK( ZLBDAS(:),MASK=GACC(:) )
    ZVEC2(:) = PACK( ZLBDAR(:),MASK=GACC(:) )
!
!        5.2.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAS)-0.00001,           &
                          XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
!
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAR)-0.00001,           &
                          XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
!
!        5.2.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
    ZVEC3(:)  = BI_LIN_INTP_V(XKER_RACCSS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGACC)
    ZZW(:)    = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
    ZVECQ5(:) = BI_LIN_INTP_V(XKER_Q_RACCSS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGACC)
    ZWQ1(:,5) = UNPACK( VECTOR=ZVECQ5(:), MASK=GACC, FIELD=0.0 )
!
!        5.2.4  raindrop accretion on the small sized aggregates:
! 		RRACCSS & QRACCSS
!
    WHERE ( GACC(:) ) 
      ZZW1(:,2) =                                            & !! coef of RRACCS
              XFRACCSS*( ZLBDAS(:)**XCXS )*ZRHOCOR(:)/(ZCOR00* ZRHODREF(:))  &
         *( XLBRACCS1/((ZLBDAS(:)**2)               ) +                  &
            XLBRACCS2/( ZLBDAS(:)    * ZLBDAR(:)    ) +                  &
            XLBRACCS3/(               (ZLBDAR(:)**2)) )/ZLBDAR(:)**4
      ZZW1(:,4) = MIN( ZRRS(:),ZZW1(:,2)*ZZW(:) )           ! RRACCSS
      ZRRS(:) = ZRRS(:) - ZZW1(:,4)
      ZRSS(:) = ZRSS(:) + ZZW1(:,4)
      ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRACCSS))
      ZWQ4(:) = XFQRACCS * ZERT(:) * ZRHOCOR(:)/(ZCOR00* ZRHODREF(:))     * &
                ZLBDAR(:)**XCXR * ZLBDAS(:)**XCXS                         * &
               (XLBQRACCS1 * ZLBDAR(:)**(-2.0 - XFR)                      + &
                XLBQRACCS2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAS(:)**(-1.0)  + &
                XLBQRACCS3 * ZLBDAR(:)**(-XFR)       * ZLBDAS(:)**(-2.0))
      ZWQ1(:,5) = ZWQ1(:,5) * ZWQ4(:)                       ! QRACCSS
    END WHERE
!
    WHERE (ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRRS(:) > ZRSMIN_ELEC(3) .AND. ZZT(:) < XTT             .AND. &
           ABS(ZQRS(:)) > XQTMIN(3) .AND. ABS(ZERT) > XERMIN) 
      ZWQ1(:,5) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,5)) ),ZQRS(:) )
      ZQRS(:) = ZQRS(:) - ZWQ1(:,5)
      ZQSS(:) = ZQSS(:) + ZWQ1(:,5)
    ENDWHERE
!
!        5.2.5  perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
    ZVEC3(:)  = BI_LIN_INTP_V(XKER_RACCS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGACC) 
    ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) 
                                                                       !! RRACCS!
!
    ZVECQ4(:) = BI_LIN_INTP_V(XKER_Q_RACCS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGACC)
    ZWQ1(:,4) = UNPACK( VECTOR=ZVECQ4(:), MASK=GACC, FIELD=0.0 )
!
!        5.2.6  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
    ZVEC3(:) = BI_LIN_INTP_V(XKER_SACCRG, IVEC2, IVEC1, ZVEC2, ZVEC1, IGACC)
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
    ZVECQ6(:) = BI_LIN_INTP_V(XKER_Q_SACCRG, IVEC2,IVEC1, ZVEC2, ZVEC1, IGACC) 
    ZWQ1(:,6) = UNPACK( VECTOR=ZVECQ6(:), MASK=GACC, FIELD=0.0 )
    ZWQ1(:,4) = ZWQ1(:,4) * ZWQ4(:)                                   ! QRACCS
!
!        5.2.7  raindrop accretion-conversion of the large sized aggregates
!               into graupeln: RRACCSG & QRACCSG and RSACCRG & QSACCRG
!
    WHERE ( GACC(:) .AND. (ZRSS(:)>0.0) ) 
      ZZW1(:,2) = MIN( ZRRS(:),ZZW1(:,2)-ZZW1(:,4) )                  ! RRACCSG
      ZZW1(:,3) = MIN( ZRSS(:),XFSACCRG*ZZW(:)*                     & ! RSACCRG
            ( ZLBDAS(:)**(XCXS-XBS) )*ZRHOCOR(:)/(ZCOR00* ZRHODREF(:))  &
           *( XLBSACCR1/((ZLBDAR(:)**2)               ) +           &
              XLBSACCR2/( ZLBDAR(:)    * ZLBDAS(:)    ) +           &
              XLBSACCR3/(               (ZLBDAS(:)**2)) )/ZLBDAR(:) )
      ZRRS(:) = ZRRS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2)+ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRACCSG))
      ZWQ1(:,4) = ZWQ1(:,4) - ZWQ1(:,5)         ! QRACCSG
      ZWQ1(:,6) = ZWQ1(:,6) * XFQRACCS * ZEST(:) *                             &
                  ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *                        &
                  ZLBDAR(:)**XCXR * ZLBDAS(:)**XCXS *                          &
                 (XLBQSACCRG1 * ZLBDAS(:)**(-2.0 - XFS) +                      &
                  XLBQSACCRG2 * ZLBDAS(:)**(-1.0 - XFS) * ZLBDAR(:)**(-1.0) +  &
                  XLBQSACCRG3 * ZLBDAS(:)**(-XFS) * ZLBDAR(:)**(-2.0)) ! QSACCR
    END WHERE
!
    WHERE (ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
           ZZT(:) < XTT             .AND. ABS(ZQGS(:)) > XQTMIN(6)) 
      ZWQ1(:,4) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,4)) ),ZQRS(:) )
      ZQRS(:) = ZQRS(:) - ZWQ1(:,4)
      ZQGS(:) = ZQGS(:) + ZWQ1(:,4)
    ENDWHERE
!
    WHERE (ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
           ZZT(:) < XTT             .AND. ABS(ZQGS(:)) > XQTMIN(6) .AND. &
           ABS(ZEST) > XESMIN) 
      ZWQ1(:,6) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,6)) ),ZQSS(:) )
      ZQSS(:) = ZQSS(:) - ZWQ1(:,6)
      ZQGS(:) = ZQGS(:) + ZWQ1(:,6)
    ENDWHERE
!
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
    DEALLOCATE( ZVECQ4 )
    DEALLOCATE( ZVECQ5 )
    DEALLOCATE( ZVECQ6 )

    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'ACC', &
                                             Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACC', &
                                             Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'ACC', &
                                             Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'ACC', &
                                             Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'ACC', &
                             Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'ACC', &
                             Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'ACC', &
                             Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
  END IF

  DEALLOCATE(GACC)
!
!*       5.3    Conversion-Melting of the aggregates: RSMLT & QSMLT
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'CMEL', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'CMEL', &
                            Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
  WHERE ((ZRST(:) > XRTMIN(5)) .AND. (ZRSS(:) > 0.0) .AND. (ZZT(:) > XTT))
    ZZW(:) = ZRVT(:) * ZPRES(:) / ((XMV / XMD) + ZRVT(:)) ! Vapor pressure
    ZZW(:) = ZKA(:) * (XTT - ZZT(:)) +                          &
            (ZDV(:) * (XLVTT + (XCPV - XCL) * (ZZT(:) - XTT)) * &
                      (XESTT - ZZW(:)) / (XRV * ZZT(:)))
!
! compute RSMLT
!
    ZZW(:) = MIN( ZRSS(:), XFSCVMG * MAX( 0.0,( -ZZW(:) *            &
                          (X0DEPS *          ZLBDAS(:)**XEX0DEPS +   &
                           X1DEPS * ZCJ(:) * ZLBDAS(:)**XEX1DEPS ) - &
                                   (ZZW1(:,1) + ZZW1(:,4)) *         &
                          (ZRHODREF(:) * XCL * (XTT - ZZT(:)))) /    &
                                            (ZRHODREF(:) * XLMTT)))
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
    ZRSS(:) = ZRSS(:) - ZZW(:)
    ZRGS(:) = ZRGS(:) + ZZW(:)
    ZWQ1(:,7) = XCOEF_RQ_S * ZQST(:) * ZZW(:) / ZRST(:)      ! QSMLT
  END WHERE
!
  WHERE (ZRST(:) > XRTMIN_ELEC(5) .AND. ZRSS(:) > ZRSMIN_ELEC(5) .AND. &
         ZRGT(:) > XRTMIN_ELEC(6) .AND. ABS(ZQST(:)) > XQTMIN(5) .AND. &
         ZZT(:) > XTT             .AND. ZRHODREF(:)*XLMTT > 0.) 
    ZWQ1(:,7) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,7)) ),ZQSS(:) )
    ZQGS(:) = ZQGS(:) + ZWQ1(:,7)
    ZQSS(:) = ZQSS(:) - ZWQ1(:,7)
  ENDWHERE

  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'CMEL', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CMEL', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'CMEL', &
                           Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'CMEL', &
                           Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  END SUBROUTINE RAIN_ICE_ELEC_FAST_RS
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_ELEC_FAST_RG
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing: RICFRRG & QICFRRG and RRCFRIG & QRCFRIG
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'CFRZ', &
                            Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'CFRZ', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'CFRZ', &
                            Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW1(:,3:4) = 0.0
  ZWQ1(:,3:4) = 0.0
  WHERE ((ZRIT(:) > XRTMIN(4)) .AND. (ZRRT(:) > XRTMIN(3)) .AND.  &
         (ZRIS(:) > 0.0)       .AND. (ZRRS(:) > 0.0))
    ZZW1(:,3) = MIN( ZRIS(:),XICFRR * ZRIT(:)                & ! RICFRRG
                                    * ZLBDAR(:)**XEXICFRR    &
                                    * ZRHOCOR(:) / ZCOR00 )
    ZZW1(:,4) = MIN( ZRRS(:),XRCFRI * ZCIT(:)                & ! RRCFRIG
                                    * ZLBDAR(:)**XEXRCFRI    &
                                    * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) )
    ZRIS(:) = ZRIS(:) - ZZW1(:,3)
    ZRRS(:) = ZRRS(:) - ZZW1(:,4)
    ZRGS(:) = ZRGS(:) + ZZW1(:,3) + ZZW1(:,4)
    ZTHS(:) = ZTHS(:) + ZZW1(:,4) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*RRCFRIG)
    ZWQ1(:,4) = XQRCFRIG * ZLBDAR(:)**XEXQRCFRIG * ZCIT(:) * &
                ZERT(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))  ! QRCFRIG
    ZWQ1(:,3) = XCOEF_RQ_I * ZQIT(:) * ZZW1(:,3) / ZRIT(:)     ! QICFRRG
  END WHERE
!
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
         ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
         ABS(ZERT) > XERMIN       .AND. ABS(ZQRT(:)) > XQTMIN(3))
    ZWQ1(:,4) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,4)) ),ZQRS(:) )
    ZQGS(:) = ZQGS(:) + ZWQ1(:,4)
    ZQRS(:) = ZQRS(:) - ZWQ1(:,4)
  ENDWHERE
!
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
         ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZRRS(:) > ZRSMIN_ELEC(3) .AND. &
         ABS(ZQIT(:)) > XQTMIN(4)) 
    ZWQ1(:,3) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,3)) ),ZQIS(:) )
    ZQGS(:) = ZQGS(:) + ZWQ1(:,3)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,3)
  ENDWHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'CFRZ', Unpack( zzw1(:,4) * ( zlsfact(:) - zlvfact(:) ) &
                                                                               * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'CFRZ', &
                                           Unpack( -zzw1(:, 4)                 * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'CFRZ', &
                                           Unpack( -zzw1(:, 3)                 * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CFRZ', &
                                           Unpack( ( zzw1(:, 3) + zzw1(:, 4) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'CFRZ', &
                           Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'CFRZ', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'CFRZ', &
                           Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
!
!*       6.2    compute the Dry growth case
!
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETG', &
                                            Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETG', &
                                            Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETG', &
                                            Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETG', &
                                            Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETG', &
                                            Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETG', &
                                            Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETG', &
                                            Unpack( zrhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'WETG', &
                            Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'WETG', &
                            Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'WETG', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'WETG', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'WETG', &
                            Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( krr == 7 ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'WETG', &
                            Unpack( zqhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW1(:,:) = 0.0
  ZWQ1(:,1:10) = 0.0
  ZWQ3(:)      = 0.0
  ZWQ4(:)      = 0.0
!
!*       6.2.1  compute RCDRYG & QCDRYG
!
  WHERE ((ZRGT(:) > XRTMIN(6)) .AND. ((ZRCT(:) > XRTMIN(2) .AND. ZRCS(:) > 0.0)))
    ZZW(:) = ZLBDAG(:)**(XCXG-XDG-2.0) * ZRHOCOR(:) / ZCOR00
    ZZW1(:,1) = MIN( ZRCS(:),XFCDRYG * ZRCT(:) * ZZW(:) )             ! RCDRYG
    ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW1(:,1) / ZRCT(:)            ! QCDRYG
  END WHERE
!
  WHERE (ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
         ABS(ZQCT(:)) > XQTMIN(2) .AND. ZRCS(:) > ZRSMIN_ELEC(2)) 
    ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
  ELSEWHERE
    ZWQ1(:,1) = 0.
  ENDWHERE
!
!*       6.2.2  compute RIDRYG & QIDRYG
!
  WHERE ((ZRGT(:) > XRTMIN(6)) .AND. ((ZRIT(:) > XRTMIN(4) .AND. ZRIS(:) > 0.0)) )
    ZZW(:) = ZLBDAG(:)**(XCXG-XDG-2.0) * ZRHOCOR(:)/ZCOR00
    ZZW1(:,2) = MIN( ZRIS(:),XFIDRYG * EXP( XCOLEXIG*(ZZT(:)-XTT) ) &
                                     * ZRIT(:) * ZZW(:) )             ! RIDRYG
   ZWQ1(:,2) = XCOEF_RQ_I * ZQIT(:) * ZZW1(:,2) / ZRIT(:)         ! QIDRYG_coal
  END WHERE
!
  WHERE (GELEC(:,2)) 
    ZWQ1(:,2) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,2)) ),ZQIS(:) ) 
  ELSEWHERE
    ZWQ1(:,2) = 0.
  ENDWHERE
!  
  CALL ELEC_IDRYG_B()                ! QIDRYG_boun
!
! Save the NI charging rate for temporal series
  XNI_IDRYG(:,:,:) = UNPACK(ZWQ1(:,3), MASK=GMICRO, FIELD=0.0)
  XNI_IDRYG(:,:,:) = XNI_IDRYG(:,:,:) * PRHODREF(:,:,:)   ! C/m3/s
!
!*       6.2.3  accretion of aggregates on the graupeln
!
  ALLOCATE(GDRY(IMICRO))
  GDRY(:) = (ZRST(:)>XRTMIN(5)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRSS(:)>0.0)
  IGDRY = COUNT( GDRY(:) )
!
  IF( IGDRY>0 ) THEN
!
!        6.2.3.1  allocations
!
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
!
    ALLOCATE( ZVECQ4(IGDRY) )
    ALLOCATE( ZVECQ5(IGDRY) )
    ALLOCATE( ZVECQ6(IGDRY) )
!
    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAUN1' .OR. &
        CNI_CHARGING == 'SAUN2' .OR. CNI_CHARGING == 'SAP98' .OR. &
        CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
        CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR' .OR. &
        CNI_CHARGING == 'GARDI') &
!
      ALLOCATE( ZAUX(IGDRY) )
!
!        6.2.3.2  select the (ZLBDAG,ZLBDAS) couplet
!
    ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
    ZVEC2(:) = PACK( ZLBDAS(:),MASK=GDRY(:) )
!
!        6.2.3.3  find the next lower indice for the ZLBDAG and for the ZLBDAS
!                 in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!                 tabulate the SDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAS)-0.00001,           &
                          XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2S ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!        6.2.3.4  perform the bilinear interpolation of the normalized
!                 SDRYG-kernel
!
! normalized SDRYG-kernel
    ZVEC3(:) = BI_LIN_INTP_V(XKER_SDRYG,   IVEC1, IVEC2, ZVEC1, ZVEC2, IGDRY)
    ZZW(:)   = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
! normalized Q-SDRYG-kernel 
    ZVECQ4(:) = BI_LIN_INTP_V(XKER_Q_SDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGDRY)
    ZWQ1(:,4) = UNPACK( VECTOR=ZVECQ4(:), MASK=GDRY, FIELD=0.0 )
!
! normalized Q-???-kernel
    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAUN1' .OR. &
        CNI_CHARGING == 'SAUN2' .OR. CNI_CHARGING == 'SAP98' .OR. &
        CNI_CHARGING == 'GARDI' .OR.                              &
        CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
        CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
      ZAUX(:)  = BI_LIN_INTP_V(XKER_Q_LIMSG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGDRY)
      ZAUX1(:) = UNPACK( VECTOR=ZAUX(:), MASK=GDRY, FIELD=0.0 )
    END IF
!
! normalized Q-SDRYG-bouncing kernel
    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'HELFA' .OR. &
        CNI_CHARGING == 'GARDI') THEN
      ZVECQ5(:)  = BI_LIN_INTP_V(XKER_Q_SDRYGB,IVEC1,IVEC2,ZVEC1,ZVEC2,IGDRY)
      ZWQ1(:,10) = UNPACK( VECTOR=ZVECQ5(:), MASK=GDRY, FIELD=0.0 )
    ELSE
      ZVECQ5(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB1,IVEC1,IVEC2,ZVEC1,ZVEC2,IGDRY)
      ZWQ3(:)   = UNPACK( VECTOR=ZVECQ5(:), MASK=GDRY, FIELD=0.0 ) ! Dvqsgmn if charge>0
      ZVECQ6(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB2,IVEC1,IVEC2,ZVEC1,ZVEC2,IGDRY)
      ZWQ4(:)   = UNPACK( VECTOR=ZVECQ6(:), MASK=GDRY, FIELD=0.0 ) ! Dvqsgmn if charge<0
    ENDIF
!
! 	6.2.3.5	 compute RSDRYG and QSDRYG = QSDRYG_coal + QSDRYG_boun
!
    WHERE( GDRY(:) )
      ZZW1(:,3) = MIN( ZRSS(:),XFSDRYG*ZZW(:)                         & ! RSDRYG
                                      * EXP( XCOLEXSG*(ZZT(:)-XTT) )  &
                    *( ZLBDAS(:)**(XCXS-XBS) )*( ZLBDAG(:)**XCXG )    &
                    * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))             &
                         *( XLBSDRYG1/( ZLBDAG(:)**2              ) + &
                            XLBSDRYG2/( ZLBDAG(:)   * ZLBDAS(:)   ) + &
                            XLBSDRYG3/(               ZLBDAS(:)**2) ) )
      ZWQ1(:,4) = ZWQ1(:,4) * XFQSDRYG *                                 &
                XCOLSG * EXP(XCOLEXSG * (ZZT(:) - XTT)) *                &
                ZEST(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *          & 
                ZLBDAG(:)**XCXG * ZLBDAS(:)**XCXS *                      &
               (XLBQSDRYG1 * ZLBDAS(:)**(-2.0-XFS) +                     &
                XLBQSDRYG2 * ZLBDAS(:)**(-1.0-XFS) * ZLBDAG(:)**(-1.0) + &
                XLBQSDRYG3 * ZLBDAS(:)**(-XFS)     * ZLBDAG(:)**(-2.0)) ! QSDRYG_coal
    END WHERE
!
    WHERE (ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ABS(ZQST(:)) > XQTMIN(5) .AND. &
           ABS(ZEST) > XESMIN)
      ZWQ1(:,4) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,4)) ),ZQSS(:) )
    ELSEWHERE
      ZWQ1(:,4) = 0.
    END WHERE
!
! QSDRYG_boun
    CALL ELEC_SDRYG_B()
!
! save the NI charging rate for temporal series
    XNI_SDRYG(:,:,:) = UNPACK(ZWQ1(:,5), MASK=GMICRO, FIELD=0.0)
    XNI_SDRYG(:,:,:) = XNI_SDRYG(:,:,:) * PRHODREF(:,:,:)  ! C/m3/
!
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
!
    DEALLOCATE( ZVECQ4 )
    DEALLOCATE( ZVECQ5 )
    DEALLOCATE( ZVECQ6 )
    IF (ALLOCATED(ZAUX)) DEALLOCATE( ZAUX )
  END IF
!
!
!*       6.2.4     accretion of raindrops on the graupeln
!
  GDRY(:) = (ZRRT(:)>XRTMIN(3)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRRS(:)>0.0)
  IGDRY = COUNT( GDRY(:) )
!
  IF( IGDRY>0 ) THEN
!
!        6.2.4.1  allocations
!
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
    ALLOCATE(ZVECQ4(IGDRY))
!
!        6.2.4.2  select the (ZLBDAG,ZLBDAR) couplet
!
    ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
    ZVEC2(:) = PACK( ZLBDAR(:),MASK=GDRY(:) )
!
!        6.2.4.3  find the next lower indice for the ZLBDAG and for the ZLBDAR
!                 in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!                 tabulate the RDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAR)-0.00001,           &
                          XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2R ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!        6.2.4.4  perform the bilinear interpolation of the normalized
!                 RDRYG-kernel
!
    ZVEC3(:) = BI_LIN_INTP_V(XKER_RDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGDRY)
    ZZW(:)   = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
    ZVECQ4(:) = BI_LIN_INTP_V(XKER_Q_RDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGDRY)
    ZWQ1(:,6) = UNPACK( VECTOR=ZVECQ4(:), MASK=GDRY, FIELD=0.0 )
!
!        6.2.4.5  compute RRDRYG and QRDRYG
!
    WHERE( GDRY(:) )
      ZZW1(:,4) = MIN( ZRRS(:),XFRDRYG*ZZW(:)                    & ! RRDRYG
                        *( ZLBDAR(:)**(-4) )*( ZLBDAG(:)**XCXG ) &
                        * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))    &
                    *( XLBRDRYG1/( ZLBDAG(:)**2              ) + &
                       XLBRDRYG2/( ZLBDAG(:)   * ZLBDAR(:)   ) + &
                       XLBRDRYG3/(               ZLBDAR(:)**2) ) )
      ZWQ1(:,6) =  ZWQ1(:,6) * XFQRDRYG *                              &
          ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *                        &
          ZERT(:) * ZLBDAG(:)**XCXG * ZLBDAR(:)**XCXR *                &
         (XLBQRDRYG1 * ZLBDAR(:)**(-2.0 - XFR) +                       &
          XLBQRDRYG2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAG(:)**(-1.0) +   &
          XLBQRDRYG3 * ZLBDAR(:)**(-XFR)       * ZLBDAG(:)**(-2.0))    ! QRDRYG
    END WHERE
!
    WHERE (ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. & 
           ZRRS(:) > ZRSMIN_ELEC(3).AND. ABS(ZERT) > XERMIN .AND.        &
           ABS(ZQRT(:)) > XQTMIN(3))
      ZWQ1(:,6) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,6)) ),ZQRS(:) )
    ELSEWHERE
      ZWQ1(:,6) = 0.
    ENDWHERE
!
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVECQ4)
  END IF
!
  ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
  DEALLOCATE(GDRY)
!
!
!*       6.3    compute the Wet growth case
!
  ZZW(:) = 0.0
  ZRWETG(:) = 0.0
  ZWQ1(:,7:9) = 0.0
!
  WHERE (ZRGT(:) > XRTMIN(6))
    ZZW1(:,5) = MIN( ZRIS(:),                                    &
                ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(ZZT(:)-XTT)) ) ) ! RIWETG
    ZZW1(:,6) = MIN( ZRSS(:),                                    &
                ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(ZZT(:)-XTT)) ) ) ! RSWETG
!
    ZZW(:) = ZRVT(:) * ZPRES(:) / ((XMV / XMD) + ZRVT(:)) ! Vapor pressure
    ZZW(:) = ZKA(:) * (XTT - ZZT(:)) +                              &
            (ZDV(:) * (XLVTT + (XCPV - XCL) * (ZZT(:) - XTT)) * &
                      (XESTT - ZZW(:)) / (XRV * ZZT(:)))
!
! compute RWETG
!
    ZRWETG(:) = MAX(0.0,                                                     &
                   (ZZW(:) * (X0DEPG *          ZLBDAG(:)**XEX0DEPG +        &
                              X1DEPG * ZCJ(:) * ZLBDAG(:)**XEX1DEPG) +       &
                   (ZZW1(:,5) + ZZW1(:,6) ) *                                &
                   (ZRHODREF(:) * (XLMTT + (XCI - XCL) * (XTT - ZZT(:))))) / &
                                (ZRHODREF(:) * (XLMTT - XCL * (XTT - ZZT(:)))))
  END WHERE
!
  WHERE (ZRGT(:) > 0.0 .AND. ZRIT(:) > 0. .AND. ZRST(:) > 0.)
    ZWQ1(:,7) = XCOEF_RQ_I * ZZW1(:,5) * ZQIT(:) / ZRIT(:)
    ZWQ1(:,8) = XCOEF_RQ_S * ZZW1(:,6) * ZQST(:) / ZRST(:)
  END WHERE
!
  WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE) 
    ZWQ1(:,7) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,7)) ),ZQIS(:) )
    ZWQ1(:,8) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,8)) ),ZQSS(:) )
  ELSEWHERE
    ZWQ1(:,7) = 0.
    ZWQ1(:,8) = 0.
  ENDWHERE
!
  WHERE (ZRGS(:) > ZRSMIN_ELEC(6) .AND. ABS(ZQRT(:)) > XQTMIN(3) .AND. &
         ZRRT(:) > XRTMIN_ELEC(3))
    ZWQ1(:,9) = XCOEF_RQ_R * ZQRT(:) * &
               (ZRWETG(:) - ZZW1(:,5) - ZZW1(:,6) - ZZW1(:,1)) / ZRRT(:)  ! QRWETG
    ZWQ1(:,9) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,9)) ),ZQRS(:) )
  ENDWHERE
!
!
!*       6.4    Select Wet or Dry case
!
  ZZW(:) = 0.0
  IF (KRR == 7) THEN
    WHERE( ZRGT(:) > XRTMIN(6)    .AND. ZZT(:) < XTT .AND.  & ! Wet
           ZRDRYG(:) >= ZRWETG(:) .AND. ZRWETG(:) > 0.0 )     ! case
      ZZW(:) = ZRWETG(:) - ZZW1(:,5) - ZZW1(:,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
      ZZW1(:,7) = MAX( 0.0,MIN( ZZW(:),ZRRS(:)+ZZW1(:,1) ) )
      ZUSW(:)   = ZZW1(:,7) / ZZW(:)
      ZZW1(:,5) = ZZW1(:,5) * ZUSW(:)
      ZZW1(:,6) = ZZW1(:,6) * ZUSW(:)
      ZRWETG(:) = ZZW1(:,7) + ZZW1(:,5) + ZZW1(:,6)
!
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRIS(:) = ZRIS(:) - ZZW1(:,5)
      ZRSS(:) = ZRSS(:) - ZZW1(:,6)
!
! assume a linear percent of conversion of graupel into hail
!
      ZRGS(:) = ZRGS(:) + ZRWETG(:)                     !     Wet growth
      ZZW(:)  = ZRGS(:) * ZRDRYG(:) / (ZRWETG(:) + ZRDRYG(:)) !        and
      ZRGS(:) = ZRGS(:) - ZZW(:)                        !   partial conversion
      ZRHS(:) = ZRHS(:) + ZZW(:)                        ! of the graupel into hail
!
      ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,7) + ZZW1(:,1) )
      ZTHS(:) = ZTHS(:) + ZZW1(:,7)*(ZLSFACT(:)-ZLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
!
      ZQCS(:) = ZQCS(:) - ZWQ1(:,1)    ! QCDRYG .equiv. QCWETG
      ZQRS(:) = ZQRS(:) - ZWQ1(:,9)
      ZQIS(:) = ZQIS(:) - ZWQ1(:,7)
      ZQSS(:) = ZQSS(:) - ZWQ1(:,8)
      ZQGS(:) = ZQGS(:) + ZWQ1(:,1) + ZWQ1(:,9) + ZWQ1(:,7) + ZWQ1(:,8)
      ZZW(:)  = ZQGS(:) * ZRDRYG(:) / (ZRWETG(:) + ZRDRYG(:)) ! partial graupel
      ZQGS(:) = ZQGS(:) - ZZW(:)                        ! charge conversion
      ZQHS(:) = ZQHS(:) + ZZW(:)                        ! into hail charge
    END WHERE
  ELSE IF( KRR == 6 ) THEN
    WHERE (ZRGT(:) > XRTMIN(6)    .AND. ZZT(:) < XTT .AND. & ! Wet
           ZRDRYG(:) >= ZRWETG(:) .AND. ZRWETG(:) > 0.0)     ! case
      ZZW(:)  = ZRWETG(:)
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRIS(:) = ZRIS(:) - ZZW1(:,5)
      ZRSS(:) = ZRSS(:) - ZZW1(:,6)
      ZRGS(:) = ZRGS(:) + ZZW(:)
!
      ZRRS(:) = ZRRS(:) - ZZW(:) + ZZW1(:,5) + ZZW1(:,6) + ZZW1(:,1)
      ZTHS(:) = ZTHS(:) + (ZZW(:)-ZZW1(:,5)-ZZW1(:,6))*(ZLSFACT(:)-ZLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
!
      ZQCS(:) = ZQCS(:) - ZWQ1(:,1)    ! QCDRYG .equiv. QCWETG
      ZQRS(:) = ZQRS(:) - ZWQ1(:,9)
      ZQIS(:) = ZQIS(:) - ZWQ1(:,7)
      ZQSS(:) = ZQSS(:) - ZWQ1(:,8)
      ZQGS(:) = ZQGS(:) + ZWQ1(:,1) + ZWQ1(:,9) + ZWQ1(:,7) + ZWQ1(:,8)
    END WHERE
  END IF

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETG', &
                                           Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETG', &
                                           Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETG', &
                                           Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETG', &
                                           Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETG', &
                                           Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETG', &
                                           Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETG', &
                                           Unpack( zrhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'WETG', &
                           Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'WETG', &
                           Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'WETG', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'WETG', &
                           Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'WETG', &
                           Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( krr == 7 ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'WETG', &
                           Unpack( zqhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DRYG', &
                                            Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DRYG', &
                                            Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'DRYG', &
                                            Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'DRYG', &
                                            Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'DRYG', &
                                            Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'DRYG', &
                                            Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'DRYG', &
                            Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'DRYG', &
                            Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'DRYG', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'DRYG', &
                            Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'DRYG', &
                            Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  WHERE (ZRGT(:) > XRTMIN(6) .AND. ZZT(:) < XTT .AND. & ! Dry
         ZRDRYG(:) < ZRWETG(:) .AND. ZRDRYG(:) > 0.0)   ! case
    ZRCS(:) = ZRCS(:) - ZZW1(:,1)
    ZRIS(:) = ZRIS(:) - ZZW1(:,2)
    ZRSS(:) = ZRSS(:) - ZZW1(:,3)
    ZRRS(:) = ZRRS(:) - ZZW1(:,4)
    ZRGS(:) = ZRGS(:) + ZRDRYG(:)
    ZTHS(:) = ZTHS(:) + (ZZW1(:,1) + ZZW1(:,4)) * (ZLSFACT(:) - ZLVFACT(:))
                                       ! f(L_f*(RCDRYG+RRDRYG))
!
    ZQCS(:) = ZQCS(:) - ZWQ1(:,1)
    ZQRS(:) = ZQRS(:) - ZWQ1(:,6)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,2) - ZWQ1(:,3)
    ZQSS(:) = ZQSS(:) - ZWQ1(:,4) - ZWQ1(:,5) 
    ZQGS(:) = ZQGS(:) + ZWQ1(:,1) + ZWQ1(:,2) + ZWQ1(:,3) + ZWQ1(:,4) & 
                                              + ZWQ1(:,5) + ZWQ1(:,6)
  END WHERE

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DRYG', &
                                           Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DRYG', &
                                           Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'DRYG', &
                                           Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'DRYG', &
                                           Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'DRYG', &
                                           Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'DRYG', &
                                           Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'DRYG', &
                           Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'DRYG', &
                           Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'DRYG', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'DRYG', &
                           Unpack( zqss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'DRYG', &
                           Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
!
! Inductive mecanism
!
  IF (LINDUCTIVE) THEN
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1), 'INCG', &
                              Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'INCG', &
                              Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if

    ZRATE_IND(:) = 0.
    GIND(:) = ZRDRYG(:) > 0. .AND. ZRDRYG(:) < ZRWETG(:) .AND. ZZT(:) < XTT
    IIND = COUNT(GIND(:))
!
    IF (IIND > 0) CALL INDUCTIVE_PROCESS
!
    XIND_RATE(:,:,:) = 0.
    XIND_RATE(:,:,:) = UNPACK(ZRATE_IND(:), MASK=GMICRO, FIELD=0.0)
    XIND_RATE(:,:,:) = XIND_RATE(:,:,:) * PRHODREF(:,:,:)  ! C/m3/s

    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1), 'INCG', &
                             Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'INCG', &
                             Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
!*       6.5    Melting of the graupeln
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2), 'GMLT', &
                            Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'GMLT', &
                            Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
  ZWQ1(:,7) = 0.0
  WHERE ((ZRGT(:) > XRTMIN(6)) .AND. (ZRGS(:) > 0.0) .AND. (ZZT(:) > XTT))
    ZZW(:) = ZRVT(:) * ZPRES(:) / ((XMV / XMD) + ZRVT(:)) ! Vapor pressure
    ZZW(:) = ZKA(:) * (XTT - ZZT(:)) +                            &
            (ZDV(:) * (XLVTT + ( XCPV - XCL ) * (ZZT(:) - XTT)) * &
                      (XESTT - ZZW(:)) / (XRV * ZZT(:)))
! compute RGMLTR
    ZZW(:)  = MIN(ZRGS(:), MAX(0.0, (-ZZW(:) *                                &
                                    (X0DEPG *          ZLBDAG(:)**XEX0DEPG +  &
                                     X1DEPG * ZCJ(:) * ZLBDAG(:)**XEX1DEPG) - &
                                    (ZZW1(:,1) + ZZW1(:,4)) *                 &
                                    (ZRHODREF(:) * XCL * (XTT - ZZT(:)))) /   &
                                    (ZRHODREF(:) * XLMTT)))
    ZRRS(:) = ZRRS(:) + ZZW(:)
    ZRGS(:) = ZRGS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) - ZZW(:) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*(-RGMLTR))
! compute QGMLTR
    ZWQ1(:,7) = XCOEF_RQ_G * ZQGT(:) * ZZW(:) / ZRGT(:)
  END WHERE
!
!
  WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND. &
         ZZT(:) > XTT             .AND. ABS(ZQGT(:)) > XQTMIN(6))
    ZWQ1(:,7) = SIGN( MIN( ABS(ZQGS(:)),ABS(ZWQ1(:,7)) ),ZQGS(:) )
    ZQRS(:) = ZQRS(:) + ZWQ1(:,7)
    ZQGS(:) = ZQGS(:) - ZWQ1(:,7)
  ENDWHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'GMLT', Unpack( -zzw(:) * ( zlsfact(:) - zlvfact(:) ) &
                                                           * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'GMLT', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'GMLT', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2), 'GMLT', &
                           Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'GMLT', &
                           Unpack( zqgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  END SUBROUTINE RAIN_ICE_ELEC_FAST_RG
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_ELEC_FAST_RH
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
  ALLOCATE( GHAIL(IMICRO) )
  GHAIL(:) = ZRHT(:) > XRTMIN(7)
  IHAIL = COUNT(GHAIL(:))
!
  IF( IHAIL>0 ) THEN
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETH', &
                                              Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETH', &
                                              Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETH', &
                                              Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETH', &
                                              Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETH', &
                                              Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETH', &
                                              Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETH', &
                                              Unpack( zrhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
!
!*       7.2    compute the Wet growth of hail
!
    WHERE (GHAIL(:))
      ZLBDAH(:) = XLBH * (ZRHODREF(:) * MAX(ZRHT(:), XRTMIN(7)))**XLBEXH
    END WHERE
!
    ZZW1(:,:) = 0.0
    WHERE (GHAIL(:) .AND. ((ZRCT(:) > XRTMIN(2) .AND. ZRCS(:) > 0.0)))
      ZZW(:) = ZLBDAH(:)**(XCXH-XDH-2.0) * ZRHOCOR(:) / ZCOR00
      ZZW1(:,1) = MIN( ZRCS(:),XFWETH * ZRCT(:) * ZZW(:) )             ! RCWETH
    END WHERE
    WHERE (GHAIL(:) .AND. ((ZRIT(:) > XRTMIN(4) .AND. ZRIS(:) > 0.0)))
      ZZW(:) = ZLBDAH(:)**(XCXH-XDH-2.0) * ZRHOCOR(:) / ZCOR00
      ZZW1(:,2) = MIN( ZRIS(:),XFWETH * ZRIT(:) * ZZW(:) )             ! RIWETH
    END WHERE
!
!*       7.2.1  accretion of aggregates on the hailstones
!
    ALLOCATE( GWET(IMICRO) )
    GWET(:) = GHAIL(:) .AND. (ZRST(:) > XRTMIN(5) .AND. ZRSS(:) > 0.0)
    IGWET = COUNT( GWET(:) )
!
    IF (IGWET > 0) THEN
!
!*       7.2.2  allocations
!
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       7.2.3  select the (ZLBDAH,ZLBDAS) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAS(:),MASK=GWET(:) )
!
!*       7.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAH)-0.00001,           &
                            XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAS)-0.00001,           &
                            XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + XWETINTP2S ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       7.2.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
      ZVEC3(:) = BI_LIN_INTP_V(XKER_SWETH, IVEC1, IVEC2, ZVEC1, ZVEC2, IGWET)
      ZZW(:)   = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,3) = MIN( ZRSS(:), XFSWETH*ZZW(:)                      & ! RSWETH
                      *( ZLBDAS(:)**(XCXS-XBS) )*( ZLBDAH(:)**XCXH )  &
                         * ZRHOCOR(:)/(ZCOR00*ZRHODREF(:))            &
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
!*       7.2.6  accretion of graupeln on the hailstones
!
    GWET(:) = GHAIL(:) .AND. (ZRGT(:)>XRTMIN(6) .AND. ZRGS(:)>0.0)
    IGWET = COUNT( GWET(:) )
!
    IF (IGWET > 0) THEN
!
!*       7.2.7  allocations
!
      ALLOCATE( ZVEC1(IGWET) )
      ALLOCATE( ZVEC2(IGWET) )
      ALLOCATE( ZVEC3(IGWET) )
      ALLOCATE( IVEC1(IGWET) )
      ALLOCATE( IVEC2(IGWET) )
!
!*       7.2.8  select the (ZLBDAH,ZLBDAG) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAG(:),MASK=GWET(:) )
!
!*       7.2.9  find the next lower indice for the ZLBDAH and for the ZLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAG)-0.00001,           &
                            XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAG)-0.00001,           &
                            XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + XWETINTP2G ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       7.2.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
      ZVEC3(:) = BI_LIN_INTP_V(XKER_GWETH, IVEC1, IVEC2, ZVEC1, ZVEC2, IGWET)
      ZZW(:)   = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE (GWET(:))
        ZZW1(:,5) = MIN( ZRGS(:),XFGWETH*ZZW(:)                       & ! RGWETH
                      *( ZLBDAG(:)**(XCXG-XBG) )*( ZLBDAH(:)**XCXH )  &
                         * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))        &
                         *( XLBGWETH1/( ZLBDAH(:)**2              ) + &
                            XLBGWETH2/( ZLBDAH(:)   * ZLBDAG(:)   ) + &
                            XLBGWETH3/(               ZLBDAG(:)**2) ) )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
    DEALLOCATE(GWET)
!
!*       7.3    compute the Wet growth of hail
!
    ZZW(:) = 0.0
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT)
      ZZW(:) = ZRVT(:) * ZPRES(:) / ((XMV / XMD) + ZRVT(:)) ! Vapor pressure
      ZZW(:) = ZKA(:) * (XTT - ZZT(:)) +                          &
              (ZDV(:) * (XLVTT + (XCPV - XCL) * (ZZT(:) - XTT)) * &
              (XESTT - ZZW(:)) / (XRV * ZZT(:)))
!
! compute RWETH
!
      ZZW(:) = MAX(0., (ZZW(:) * (X0DEPH *          ZLBDAH(:)**XEX0DEPH + &
                          X1DEPH * ZCJ(:) * ZLBDAH(:)**XEX1DEPH) +        &
               (ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5) ) *                     &
               (ZRHODREF(:) * (XLMTT + (XCI - XCL) * (XTT - ZZT(:))))) /  &
                    (ZRHODREF(:) * (XLMTT - XCL * (XTT - ZZT(:)))))
!
      ZZW1(:,6) = MAX( ZZW(:) - ZZW1(:,2) - ZZW1(:,3) - ZZW1(:,5), 0. ) ! RCWETH+RRWETH
    END WHERE
!
    ZUSW(:)   = 0.
!
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZZW1(:,6) /= 0.0)
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
      ZZW1(:,4) = MAX( 0.0,MIN( ZZW1(:,6),ZRRS(:)+ZZW1(:,1) ) )
      ZUSW(:)   = ZZW1(:,4) / ZZW1(:,6)
      ZZW1(:,2) = ZZW1(:,2)*ZUSW(:)
      ZZW1(:,3) = ZZW1(:,3)*ZUSW(:)
      ZZW1(:,5) = ZZW1(:,5)*ZUSW(:)
      ZZW(:)    = ZZW1(:,4) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5)
!
!*       7.1.6  integrate the Wet growth of hail
!
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRIS(:) = ZRIS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) - ZZW1(:,5)
      ZRHS(:) = ZRHS(:) + ZZW(:)
      ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,4) + ZZW1(:,1) )
      ZRRS(:) = ZRRS(:) - ZZW1(:,4)
      ZTHS(:) = ZTHS(:) + (ZZW1(:,4)+ZZW1(:,1))*(ZLSFACT(:)-ZLVFACT(:))
                                 ! f(L_f*(RCWETH+RRWETH))
    END WHERE
!
    ZWQ1(:,:) = 0.0
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZRCT(:) > XRTMIN_ELEC(2))
      ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW1(:,1) / ZRCT(:)
      ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
    END WHERE
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZRIT(:) > XRTMIN_ELEC(4))
      ZWQ1(:,2) = XCOEF_RQ_I * ZQIT(:) * ZZW1(:,2) / ZRIT(:)
      ZWQ1(:,2) = SIGN( MIN( ABS(ZQIS(:)),ABS(ZWQ1(:,2)) ),ZQIS(:) )
    END WHERE
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZRST(:) > XRTMIN_ELEC(5))
      ZWQ1(:,3) = XCOEF_RQ_S * ZQST(:) * ZZW1(:,3) / ZRST(:)
      ZWQ1(:,3) = SIGN( MIN( ABS(ZQSS(:)),ABS(ZWQ1(:,3)) ),ZQSS(:) )
    END WHERE
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZRGT(:) > XRTMIN_ELEC(6))
      ZWQ1(:,5) = XCOEF_RQ_G * ZQGT(:) * ZZW1(:,5) / ZRGT(:)
      ZWQ1(:,5) = SIGN( MIN( ABS(ZQGS(:)),ABS(ZWQ1(:,5)) ),ZQGS(:) )
    END WHERE
    WHERE (GHAIL(:) .AND. ZZT(:) < XTT .AND. ZRRT(:) > XRTMIN_ELEC(3))
      ZWQ1(:,4) = XCOEF_RQ_R * ZQRT(:) * ZZW1(:,4) / ZRRT(:)
      ZWQ1(:,4) = SIGN( MIN( ABS(ZQRS(:)),ABS(ZWQ1(:,4)) ),ZQRS(:) )
    END WHERE
!
    ZQCS(:) = ZQCS(:) - ZWQ1(:,1)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,2)
    ZQSS(:) = ZQSS(:) - ZWQ1(:,3)
    ZQGS(:) = ZQGS(:) - ZWQ1(:,5)
    ZQRS(:) = ZQRS(:) - ZWQ1(:,4)
    ZQHS(:) = ZQHS(:) + ZWQ1(:,1) + ZWQ1(:,2) + ZWQ1(:,3) + ZWQ1(:,4) + ZWQ1(:,5)

    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETH', &
                                             Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETH', &
                                             Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETH', &
                                             Unpack( zrrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETH', &
                                             Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETH', &
                                             Unpack( zrss(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETH', &
                                             Unpack( zrgs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETH', &
                                             Unpack( zrhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'WETH', &
                             Unpack( -zwq1(:, 1) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'WETH', &
                             Unpack( -zwq1(:, 4) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'WETH', &
                             Unpack( -zwq1(:, 2) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 4 ), 'WETH', &
                             Unpack( -zwq1(:, 3) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 5 ), 'WETH', &
                             Unpack( -zwq1(:, 5) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'WETH',                      &
                             Unpack( ( zwq1(:, 1) + zwq1(:, 2) + zwq1(:, 3) + zwq1(:, 4) + zwq1(:, 5) ) &
                                       * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
  IF (IHAIL > 0) THEN
!
!*       7.5    Melting of the hailstones
!
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'HMLT', &
                              Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'HMLT', &
                              Unpack( zqhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if

    ZZW(:) = 0.0
    ZWQ1(:,7) = 0.
!
    WHERE (GHAIL(:) .AND. (ZRHS(:) > 0.0) .AND. (ZZT(:) > XTT))
      ZZW(:) = ZRVT(:) * ZPRES(:) / ((XMV / XMD) + ZRVT(:)) ! Vapor pressure
      ZZW(:) = ZKA(:) * (XTT - ZZT(:)) +                          &
             ( ZDV(:) * (XLVTT + (XCPV - XCL) * (ZZT(:) - XTT)) * &
                        (XESTT - ZZW(:)) / (XRV * ZZT(:)))
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
! compute QHMLTR
      ZWQ1(:,7) = XCOEF_RQ_H * ZQHT(:) * ZZW(:) / ZRHT(:)
      ZTHS(:) = ZTHS(:) - ZZW(:) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*(-RHMLTR))
    END WHERE
!
    WHERE (ZRHT(:) > XRTMIN_ELEC(7) .AND. ZRHS(:) > ZRSMIN_ELEC(7) .AND. &
           ZZT(:) > XTT             .AND. ABS(ZQHT(:)) > XQTMIN(7))
      ZWQ1(:,7) = SIGN( MIN( ABS(ZQHS(:)),ABS(ZWQ1(:,7)) ),ZQHS(:) )
      ZQRS(:) = ZQRS(:) + ZWQ1(:,7)
      ZQHS(:) = ZQHS(:) - ZWQ1(:,7)
    END WHERE

    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HMLT', &
               Unpack( -zzw(:) * ( zlsfact(:) - zlvfact(:) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'HMLT', &
                                             Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'HMLT', &
                                             Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 2 ), 'HMLT', &
                             Unpack( zqrs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 6 ), 'HMLT', &
                             Unpack( zqhs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
  DEALLOCATE(GHAIL)
!
  END SUBROUTINE RAIN_ICE_ELEC_FAST_RH
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_ELEC_FAST_RI
!
!*       0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       7.1    cloud ice melting
!
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'IMLT', &
                                            Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'IMLT', &
                                            Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'IMLT', &
                                            Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'IMLT', &
                            Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'IMLT', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
  ZWQ1(:,1) = 0.0
  WHERE ((ZRIS(:) > 0.0) .AND. (ZZT(:) > XTT))
    ZZW(:)  = ZRIS(:)
    ZRCS(:) = ZRCS(:) + ZRIS(:)
    ZTHS(:) = ZTHS(:) - ZRIS(:) * (ZLSFACT(:) - ZLVFACT(:)) ! f(L_f*(-RIMLTC))
    ZRIS(:) = 0.0
    ZCIT(:) = 0.0
    ZQCS(:) = ZQCS(:) + ZQIS(:)
    ZQIS(:) = 0.
  END WHERE

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'IMLT', &
                                           Unpack( zths(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'IMLT', &
                                           Unpack( zrcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'IMLT', &
                                           Unpack( zris(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'IMLT', &
                           Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'IMLT', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'BERFI', &
                            Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'BERFI', &
                            Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  ZZW(:) = 0.0
  ZWQ1(:,1)  = 0.0
  WHERE ((ZRCS(:) > 0.0) .AND. (ZSSI(:) > 0.0) .AND.      &
         (ZRIT(:) > XRTMIN(4)) .AND. (ZCIT(:) > 0.0) .AND. &
          ZRCT(:) > 0.)
    ZZW(:) = MIN(1.E8,XLBI*( ZRHODREF(:)*ZRIT(:)/ZCIT(:) )**XLBEXI) ! Lbda_i
    ZZW(:) = MIN( ZRCS(:),( ZSSI(:) / (ZRHODREF(:)*ZAI(:)) ) * ZCIT(:) * &
                  ( X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0) ) )
    ZRCS(:) = ZRCS(:) - ZZW(:)
    ZRIS(:) = ZRIS(:) + ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCBERI))
!
    ZWQ1(:,1) = XCOEF_RQ_C * ZQCT(:) * ZZW(:) / ZRCT(:)
  END WHERE
!
  WHERE (ZRCS(:) > 0.0            .AND. ZSSI(:) > 0.0 .AND.  &
         ZRIT(:) > 0.0            .AND. ZCIT(:) > 0.0 .AND.  &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ABS(ZQCT(:)) > XQTMIN(2))
    ZWQ1(:,1) = SIGN( MIN( ABS(ZQCS(:)),ABS(ZWQ1(:,1)) ),ZQCS(:) )
    ZQIS(:) = ZQIS(:) + ZWQ1(:,1) 
    ZQCS(:) = ZQCS(:) - ZWQ1(:,1)
  ENDWHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'BERFI', &
              Unpack( zzw(:) * ( zlsfact(:) - zlvfact(:) ) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'BERFI', &
                                           Unpack( -zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'BERFI', &
                                           Unpack(  zzw(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  if ( lbudget_sv ) then
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'BERFI', &
                           Unpack( zqcs(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'BERFI', &
                           Unpack( zqis(:) * zrhodj(:), mask = gmicro(:, :, :), field = 0. ) )
  end if

  END SUBROUTINE RAIN_ICE_ELEC_FAST_RI
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE COMPUTE_LBDA(ZRR, ZRS, ZRG, ZRH)
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN) :: ZRR, ZRS, ZRG 
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: ZRH
!
!
!*      1.      COMPUTE LAMBDA
!               --------------
!
  ZLBDAR(:) = 0.0 
  ZLBDAS(:) = 0.0
  ZLBDAG(:) = 0.0
!
  WHERE( ZRR(:) > 0.0 )
    ZLBDAR(:) = XLBR * (ZRHODREF(:) * MAX(ZRR(:), XRTMIN(3)))**XLBEXR
  END WHERE
!
  WHERE ( ZRS(:) > 0.0 )
    ZLBDAS(:) = MIN( XLBDAS_MAX,                                           &
                     XLBS * (ZRHODREF(:) * MAX(ZRS(:), XRTMIN(5)))**XLBEXS )
  END WHERE
!
  WHERE ( ZRG(:) > 0.0 )
    ZLBDAG(:) = XLBG * (ZRHODREF(:) * MAX( ZRG(:), XRTMIN(6)))**XLBEXG
  END WHERE
!
  IF (PRESENT(ZRH)) THEN
    ZLBDAH(:) = 0.0
    WHERE ( ZRH(:) > 0.0 )
      ZLBDAH(:) = XLBH * (ZRHODREF(:) * MAX( ZRH(:), XRTMIN(7)))**XLBEXH
    END WHERE
  END IF
!
END SUBROUTINE COMPUTE_LBDA
!
!------------------------------------------------------------------------------
!
SUBROUTINE ELEC_UPDATE_QD(ZDUM, ZER, ZEI, ZES, ZEG, ZQR, ZQI, ZQS, ZQG,  &
                          ZRR, ZRI, ZRS, ZRG,                      &
                          ZEH, ZQH, ZRH, ZEC, ZQC, ZRC)
!
!	Purpose : update the parameter e_x in the relation q_x = e_x d**f_x
!                 e_x = q_x/(N_x *  M(f_x))
!
!
!*	0.	DECLARATIONS
!          	------------
!
IMPLICIT NONE
!
REAL, INTENT(IN)                :: ZDUM               ! =1. if mixing ratio
                                                      ! =timestep if source
REAL, DIMENSION(:), INTENT(IN)  :: ZQR, ZQI, ZQS, ZQG ! V. C.
REAL, DIMENSION(:), INTENT(IN)  :: ZRR, ZRI, ZRS, ZRG ! mixing ratio
REAL, DIMENSION(:), INTENT(OUT) :: ZER, ZEI, ZES, ZEG ! Coef of the charge diameter relation
REAL, DIMENSION(:), INTENT(IN),  OPTIONAL :: ZQH   ! hail
REAL, DIMENSION(:), INTENT(IN),  OPTIONAL :: ZRH   ! hail
REAL, DIMENSION(:), INTENT(OUT), OPTIONAL :: ZEH   ! hail
!
REAL, DIMENSION(:), INTENT(IN),  OPTIONAL :: ZQC   ! V. C. for droplets
REAL, DIMENSION(:), INTENT(IN),  OPTIONAL :: ZRC   ! mixing ration for droplets
REAL, DIMENSION(:), INTENT(OUT), OPTIONAL :: ZEC   ! Coef of the charge diameter relation for droplets
REAL, DIMENSION(SIZE(XRTMIN)) :: ZRTMIN_E
!
!
!*       1.     UPDATE E_x
!               ---------- 
!
IF (PRESENT(ZEC))  ZEC(:) = 0.
ZER(:) = 0.
ZEI(:) = 0.
ZES(:) = 0.
ZEG(:) = 0.
ZRTMIN_E(:) = XRTMIN(:) / ZDUM
!
!*       1.1   for cloud droplets
!
IF (PRESENT(ZEC) .AND. PRESENT(ZQC) .AND. PRESENT(ZRC)) THEN 
  WHERE (ZRC(:) > ZRTMIN_E(2)) 
    ZEC(:) = ZDUM * ZRHODREF(:) * ZQC(:) / XFQUPDC
    ZEC(:) = SIGN( MIN(ABS(ZEC(:)), XECMAX), ZEC(:))
  ENDWHERE
END IF
!
!*       1.2   for raindrops
!
WHERE (ZRR(:) > ZRTMIN_E(3) .AND. ZLBDAR(:) > 0.)
  ZER(:) = ZDUM * ZRHODREF(:) * ZQR(:) / (XFQUPDR * ZLBDAR(:)**(XCXR - XFR))
  ZER(:) = SIGN( MIN(ABS(ZER(:)), XERMAX), ZER(:))
ENDWHERE
!
!*       1.3   for ice crystals
!
WHERE (ZRI(:) > ZRTMIN_E(4) .AND. ZCIT(:) > 0.0) 
  ZEI(:) = ZDUM * ZRHODREF(:) * ZQI(:) /                       &
           ((ZCIT**(1 - XEXFQUPDI)) * XFQUPDI * (ZRHODREF(:) * &
           ZDUM * ZRI(:))**XEXFQUPDI)
  ZEI(:) = SIGN( MIN(ABS(ZEI(:)), XEIMAX), ZEI(:))
ENDWHERE
!
!*       1.4   for snow
!
WHERE (ZRS(:) > ZRTMIN_E(5) .AND. ZLBDAS(:) > 0.) 
  ZES(:) = ZDUM * ZRHODREF(:) * ZQS(:) / (XFQUPDS * ZLBDAS(:)**(XCXS - XFS))
  ZES(:) = SIGN( MIN(ABS(ZES(:)), XESMAX), ZES(:))
ENDWHERE
!
!*       1.5   for graupel
!
WHERE (ZRG(:) > ZRTMIN_E(6).AND. ZLBDAG(:) > 0.)  
  ZEG(:) = ZDUM * ZRHODREF(:) * ZQG(:) / (XFQUPDG * ZLBDAG(:)**(XCXG - XFG))
  ZEG(:) = SIGN( MIN(ABS(ZEG(:)), XEGMAX), ZEG(:))
ENDWHERE
!
!*       1.6   for hail
!
IF (PRESENT(ZEH) .AND. PRESENT(ZQH) .AND. PRESENT(ZRH)) THEN 
  ZEH(:) = 0.
  WHERE (ZRH(:) > ZRTMIN_E(7).AND. ZLBDAH(:) > 0.)  
    ZEH(:) = ZDUM * ZRHODREF(:) * ZQH(:) / (XFQUPDH * ZLBDAH(:)**(XCXH - XFH))
    ZEH(:) = SIGN( MIN(ABS(ZEH(:)), XEHMAX), ZEH(:))
  ENDWHERE
END IF
!
END  SUBROUTINE ELEC_UPDATE_QD
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_INI_NI_PROCESS
!
!	Purpose : initialization for the non-inductive charging process
!
!  GELEC(:,1) : logical variable for Ice-Snow process --> ELEC_IAGGS_B
!                               from RAIN_ICE_ELEC_SLOW routine
!  GELEC(:,2) : logical variable for Ice-Graupel process --> ELEC_IDRYG_B
!                               from RAIN_ICE_ELEC_FAST_RG
!  GELEC(:,3) : logical variable for Snow-Graupel process --> ELEC_SDRYG_B
!                               from RAIN_ICE_ELEC_FAST_RG
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
!
!*      1.      Gardiner et al. (1985)
!               ----------------------
!
  IF (CNI_CHARGING == 'GARDI') THEN
    ZDELTALWC(:) = 0.
    ZFT(:) = 0.
!
    GELEC(:,3) = ZZT(:) > (XTT - 40.) .AND. ZZT(:) < XTT
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
    WHERE (GELEC(:,4))
      ZFT(:) = - 1.7E-5 * ((-21 / (XQTC - XTT)) * (ZZT(:) - XTT))**3   &
               - 0.003  * ((-21 / (XQTC - XTT)) * (ZZT(:) - XTT))**2   &
               - 0.05   * ((-21 / (XQTC - XTT)) * (ZZT(:) - XTT))      &
               + 0.13
!
      ZDELTALWC(:) = (ZRCT(:) * ZRHODREF(:) * 1.E3) - XLWCC   ! (g m^-3)
    ENDWHERE
  ENDIF
!
!
!*      2.      Saunders et al. (1991)
!               ----------------------
!
!*      2.1     common to SAUN1 and SAUN2
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
    ZDQLWC(:) = 0.
    ZEW(:)    = 0.
!
! positive case is the default value
    ZFQIAGGS(:)   = XFQIAGGSP
    ZFQIDRYGBS(:) = XFQIDRYGBSP
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
    ZSAUNIM(:)  = XIMP     !3.76
    ZSAUNIN(:)  = XINP     !2.5
    ZSAUNSK(:)  = XSKP     !52.8
    ZSAUNSM(:)  = XSMP     !0.44
    ZSAUNSN(:)  = XSNP     !2.5
!
! LWC_crit
    ZLWCC(:) = MIN( MAX( -0.49 + 6.64E-2*(XTT-ZZT(:)),0.22 ),1.1 )   ! (g m^-3)
!
! Mansell et al. (2005, JGR): droplet collection efficiency of the graupel ~ 0.6-1.0
    ZEW(:)   = 0.8 * ZRCT(:) * ZRHODREF(:) * 1.E3   ! (g m^-3)
!
    GELEC(:,3) = ZZT(:) > (XTT - 40.) .AND. ZZT(:) <= XTT .AND. &
                 ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
    ALLOCATE (GSAUN(IMICRO))
    GSAUN(:) = .FALSE.
!
! For temperature lower than -30C and higher than -40C, value of q at -30C
    GSAUN(:) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
    IGSAUN = COUNT (GSAUN(:))
!
    IF (IGSAUN > 0) THEN
      CALL ELEC_INI_NI_SAUNQ(ZEW, ZDQLWC)
!
      WHERE (ZDQLWC(:) < 0.)
        ZFQIAGGS(:)   = XFQIAGGSN
        ZFQIDRYGBS(:) = XFQIDRYGBSN
        ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
        ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
        ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
        ZSAUNIM(:) = XIMN      !2.54
        ZSAUNIN(:) = XINN      !2.8
        ZSAUNSK(:) = XSKN      !24.
        ZSAUNSM(:) = XSMN      !0.5
        ZSAUNSN(:) = XSNN      !2.8
      ENDWHERE
    ENDIF
!
    DEALLOCATE( GSAUN )
  END IF
!
!
!*      3.      Saunders and Peck (1998)
!
  IF (CNI_CHARGING == 'SAP98') THEN
    ZRAR_CRIT(:) = 0.
!
! compute the critical rime accretion rate
    WHERE (ZZT(:) <= XTT .AND. ZZT(:) >= (XTT - 23.7)) ! Original from SAP98
      ZRAR_CRIT(:) = 1.0 + 7.93E-2 * (ZZT(:) - XTT) +    &
                           4.48E-2 * (ZZT(:) - XTT)**2 + &
                           7.48E-3 * (ZZT(:) - XTT)**3 + &
                           5.47E-4 * (ZZT(:) - XTT)**4 + &
                           1.67E-5 * (ZZT(:) - XTT)**5 + &
                           1.76E-7 * (ZZT(:) - XTT)**6
    END WHERE
!
    WHERE (ZZT(:) < (XTT - 23.7) .AND. ZZT(:) > (XTT - 40.)) ! Added by Mansell
      ZRAR_CRIT(:) = 3.4 * (1.0 - (ABS(ZZT(:) - XTT + 23.7) / & ! et al. (2005)
                     (-23.7 + 40.))**3.)
    END WHERE
!
    GELEC(:,3) = ZZT(:) >= (XTT - 40.) .AND. ZZT(:) <= XTT
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
!+++++++++ I - G collisions +++++++++
    ZSAUNIM_IG(:) = 0.
!
! positive case is the default value
    ZSAUNIM_IG(:) = XIMP
    ZSAUNIN_IG(:) = XINP
!
! Compute the Rime Accretion Rate 
    ZRAR(:) = 0.
    ZVGMEAN(:) = 0.
    WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
      ZVGMEAN(:) = XVGCOEF * ZRHOCOR(:) * ZLBDAG(:)**(-XDG)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVGMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.1 
    GELEC(:,4) = GELEC(:,2)
!
    IF (COUNT(GELEC(:,4)) .GT. 0) THEN
!
! compute the coefficients for I-G collisions
      CALL ELEC_INI_NI_SAP98 (ZRAR, ZDQRAR_IG)
!
      WHERE (ZDQRAR_IG(:) < 0.)
        ZSAUNIM_IG(:) = XIMN
        ZSAUNIN_IG(:) = XINN
      ENDWHERE
    ENDIF
!
!+++++++++ I - S collisions +++++++++
    ZDQRAR_IS(:)  = 0.
!
! positive case is the default value
    ZSAUNIM_IS(:) = XIMP
    ZSAUNIN_IS(:) = XINP
!
! Compute the Rime Accretion Rate 
    ZRAR(:) = 0.
    ZVSMEAN(:) = 0.
!
    WHERE (ZLBDAS(:) > 0. .AND. ZRCT(:) > 0.)
      ZVSMEAN(:) = XVSCOEF * ZRHOCOR(:) * ZLBDAS(:)**(-XDS)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVSMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.1
    GELEC(:,4) = GELEC(:,1)
!
    IF (COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for I-S collisions
      CALL ELEC_INI_NI_SAP98 (ZRAR, ZDQRAR_IS)
!
      WHERE (ZDQRAR_IS(:) < 0.)
        ZSAUNIM_IS(:) = XIMN
        ZSAUNIN_IS(:) = XINN
      ENDWHERE
    ENDIF
!
!+++++++++ S - G collisions +++++++++
    ZDQRAR_SG(:)  = 0.
!
! positive case is the default value
    ZSAUNSK_SG(:) = XSKP
    ZSAUNSM_SG(:) = XSMP
    ZSAUNSN_SG(:) = XSNP
!
! Compute the Rime Accretion Rate 
    ZRAR(:) = 0.
!
    WHERE (ZVSMEAN(:) > 0. .AND. ZVGMEAN(:) > 0.)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ABS(ZVGMEAN(:) - ZVSMEAN(:)) * 1.E3
    END WHERE
!
    GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.1 
    GELEC(:,4) = GELEC(:,3)
!
    IF( COUNT(GELEC(:,4)) .GT. 0) THEN
!
! compute the coefficients for S-G collisions
      CALL ELEC_INI_NI_SAP98 (ZRAR, ZDQRAR_SG)
!
      WHERE (ZDQRAR_SG(:) < 0.)
        ZSAUNSK_SG(:) = XSKN
        ZSAUNSM_SG(:) = XSMN
        ZSAUNSN_SG(:) = XSNN
      ENDWHERE
    ENDIF
  END IF
!
!*      4.      Brooks et al. (1997) without / with anomalies 
!
  IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN

    ALLOCATE (GSAUN(IMICRO))
!
! compute the critical rime accretion rate
    WHERE (ZZT(:) > (XTT - 10.7))
      ZRAR_CRIT(:) = 0.66
    END WHERE
    WHERE (ZZT(:) <= (XTT - 10.7) .AND. ZZT(:) >= (XTT - 23.7))
      ZRAR_CRIT(:) = -1.47 - 0.2 * (ZZT(:) - XTT)
    END WHERE
    WHERE (ZZT(:) < (XTT - 23.7) .AND. ZZT(:) > (XTT - 40.))
      ZRAR_CRIT(:) = 3.3
    END WHERE
!
    GELEC(:,3) = ZZT(:) > (XTT - 40.) .AND. ZZT(:) <= XTT .AND. &
                 ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
!+++++++++ I - S collisions +++++++++
    ZDQRAR_IS(:)  = 0.
!
! positive case is the default value
    ZSAUNIM_IS(:) = XIMP
    ZSAUNIN_IS(:) = XINP
!
    GSAUN(:) = .FALSE.
!
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
    ZVSMEAN(:) = 0.
!
    WHERE (ZLBDAS(:) > 0. .AND. ZRCT(:) > 0.)
      ZVSMEAN(:) = XVSCOEF * ZRHOCOR(:) * ZLBDAS(:)**(-XDS)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVSMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.1 
    GSAUN(:) = GELEC(:,1)
    IGSAUN = COUNT (GSAUN(:))
!
    IF (IGSAUN .GT. 0) THEN
      ZEW(:) = ZRAR(:) / 3.
!
      CALL ELEC_INI_NI_SAUNQ (ZEW, ZDQRAR_IS)
!
      WHERE (ZDQRAR_IS(:) < 0.)
        ZSAUNIM_IS(:) = XIMN
        ZSAUNIN_IS(:) = XINN
      ENDWHERE
    ENDIF
!
!+++++++++ I - G collisions +++++++++
    ZDQRAR_IG(:)  = 0.
!
! positive case is the default value
    ZSAUNIM_IG(:) = XIMP
    ZSAUNIN_IG(:) = XINP
!
    GSAUN(:) = .FALSE.
!
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
    ZVGMEAN(:) = 0.
!
    WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
      ZVGMEAN(:) = XVGCOEF * ZRHOCOR(:) * ZLBDAG(:)**(-XDG)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVGMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.1 
    GSAUN(:) = GELEC(:,2)
    IGSAUN = COUNT (GSAUN(:))
!
    IF (IGSAUN .GT. 0) THEN
      ZEW(:) = ZRAR(:) / 3.
      CALL ELEC_INI_NI_SAUNQ (ZEW, ZDQRAR_IG)
!
      WHERE (ZDQRAR_IG(:) < 0.)
        ZSAUNIM_IG(:) = XIMN
        ZSAUNIN_IG(:) = XINN
      ENDWHERE
    ENDIF
!
!+++++++++ S - G collisions +++++++++
    ZDQRAR_SG(:)  = 0.
!
! positive case is the default value
    ZSAUNSK_SG(:) = XSKP
    ZSAUNSM_SG(:) = XSMP
    ZSAUNSN_SG(:) = XSNP
!
    GSAUN(:) = .FALSE.
!
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
!
    WHERE (ZVSMEAN(:) > 0. .AND. ZVGMEAN(:) > 0.)
      ZRAR(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * ABS(ZVGMEAN(:) - ZVSMEAN(:)) * 1.E3
    END WHERE
!
    GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.1 
    GSAUN(:) = GELEC(:,3)
    IGSAUN = COUNT (GSAUN(:))
!
    IF (IGSAUN .GT. 0) THEN
      ZEW(:) = ZRAR(:) / 3.
      CALL ELEC_INI_NI_SAUNQ (ZEW, ZDQRAR_SG)
!
      WHERE (ZDQRAR_SG(:) < 0.)
        ZSAUNSK_SG(:) = XSKN
        ZSAUNSM_SG(:) = XSMN
        ZSAUNSN_SG(:) = XSNN
      ENDWHERE
    ENDIF
!
    DEALLOCATE( GSAUN )
  END IF
!
!
!*      5.      Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    ZDQLWC(:) = 0.
!
    ZEW(:) = ZRCT(:) * ZRHODREF(:) * 1.E3      ! (g m^-3)
!
    GELEC(:,3) = ZZT(:) > (XTT - 40.) .AND. ZZT(:) <= XTT .AND. &
                 ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
    ALLOCATE (GTAKA(IMICRO))
    GTAKA(:) = .FALSE.
!
! For temperature lower than -30C and higher than -40C, value of q at -30C
    GTAKA(:) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
    IGTAKA = COUNT (GTAKA(:))
!
    IF (IGTAKA > 0) THEN
      CALL ELEC_INI_NI_TAKAH(ZEW, ZDQLWC, XMANSELL)
    ENDIF
!
    DEALLOCATE( GTAKA )
  ENDIF
!
!
!*      6.    Takahashi with EW (Tsenova and Mitzeva, 2009)
!
  IF (CNI_CHARGING == 'TEEWC') THEN
    ZDQLWC(:) = 0.
!
! positive case is the default value
    ZFQIAGGS(:)   = XFQIAGGSP_TAK
    ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
    ZSAUNIM(:)  = XIMP     !3.76
    ZSAUNIN(:)  = XINP     !2.5
    ZSAUNSK(:)  = XSKP_TAK !6.5
    ZSAUNSM(:)  = XSMP     !0.44
    ZSAUNSN(:)  = XSNP     !2.5
!
! Compute the effective water content
    ZEW(:) = 0.
!
    WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
      ZEW(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * 1.E3 
    END WHERE
!
    GELEC(:,3) = ZZT(:) >= (XTT - 40.) .AND. ZZT(:) <= XTT .AND.  &
                 ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
    ALLOCATE (GTAKA(IMICRO))
    GTAKA(:) = .FALSE.
!
! For temperature lower than -30C and higher than -40C, value of q at -30C
    GTAKA(:) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
    IGTAKA = COUNT (GTAKA(:))
!
    IF (IGTAKA > 0) THEN
      CALL ELEC_INI_NI_TAKAH(ZEW, ZDQLWC, XTAKA_TM)
!
      WHERE (ZDQLWC(:) < 0.)
        ZFQIAGGS(:)   = XFQIAGGSN_TAK
        ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
        ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
        ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
        ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
        ZSAUNIM(:) = XIMN      !2.54
        ZSAUNIN(:) = XINN      !2.8
        ZSAUNSK(:) = XSKN_TAK  !2.0
        ZSAUNSM(:) = XSMN      !0.5
        ZSAUNSN(:) = XSNN      !2.8
      ENDWHERE
    ENDIF
!
    DEALLOCATE( GTAKA )
  ENDIF
!
!
!*      7.    Takahashi with RAR (Tsenova and Mitzeva, 2011)
!
  IF (CNI_CHARGING == 'TERAR') THEN
!
    ALLOCATE (GTAKA(IMICRO))
!
    GELEC(:,3) = ZZT(:) >= (XTT - 40.) .AND. ZZT(:) <= XTT
    GELEC(:,1) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
    GELEC(:,2) = GELEC(:,3) .AND.                                              &
                 ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
    GELEC(:,3) = GELEC(:,3) .AND.                                              &
                 ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
                 ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
                 ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
!+++++++++ I - S collisions +++++++++
    ZDQRAR_IS(:) = 0.
!
! positive case is the default value
    ZSAUNIM_IS(:) = XIMP
    ZSAUNIN_IS(:) = XINP
!
    GTAKA(:) = .FALSE.
! 
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
    ZVSMEAN(:) = 0.
!
    WHERE (ZLBDAS(:) > 0. .AND. ZRCT(:) > 0.)
      ZVSMEAN(:) = XVSCOEF * ZRHOCOR(:) * ZLBDAS(:)**(-XDS)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVSMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
    GTAKA(:) = GELEC(:,1)
!
    IGTAKA = COUNT (GTAKA(:))
!
    IF (IGTAKA > 0) THEN
      ZEW(:) = ZRAR(:) / 8.
      CALL ELEC_INI_NI_TAKAH(ZEW, ZDQRAR_IS, XTAKA_TM)
!
      WHERE (ZDQRAR_IS(:) < 0.)
        ZSAUNIM_IS(:) = XIMN
        ZSAUNIN_IS(:) = XINN
      ENDWHERE
    END IF
!
!
!+++++++++ I - G collisions +++++++++
    ZDQRAR_IG(:) = 0.
!
! positive case is the default value
    ZSAUNIM_IG(:) = XIMP
    ZSAUNIN_IG(:) = XINP
!
    GTAKA(:) = .FALSE.
! 
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
    ZVGMEAN(:) = 0.
!
    WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
      ZVGMEAN(:) = XVGCOEF * ZRHOCOR(:) * ZLBDAG(:)**(-XDG)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVGMEAN(:) * 1.E3
    END WHERE
!
    GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
    GTAKA(:) = GELEC(:,2)
!
    IGTAKA = COUNT (GTAKA(:))
!
    IF (IGTAKA > 0) THEN
      ZEW(:) = ZRAR(:) / 8.
      CALL ELEC_INI_NI_TAKAH(ZEW, ZDQRAR_IG, XTAKA_TM)
!
      WHERE (ZDQRAR_IG(:) < 0.)
        ZSAUNIM_IG(:) = XIMN
        ZSAUNIN_IG(:) = XINN
      ENDWHERE
    ENDIF
!
!+++++++++ S - G collisions +++++++++
    ZDQRAR_SG(:) = 0.
!
! positive case is the default value
    ZSAUNSK_SG(:) = XSKP_TAK
    ZSAUNSM_SG(:) = XSMP
    ZSAUNSN_SG(:) = XSNP
!
    GTAKA(:) = .FALSE.
!
! Compute the Rime Accretion Rate
    ZRAR(:) = 0.
!
    WHERE (ZVSMEAN(:) > 0. .AND. ZVGMEAN(:) > 0.)
      ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ABS(ZVGMEAN(:) - ZVSMEAN(:)) * 1.E3
    END WHERE
!
    GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80
    GTAKA(:) = GELEC(:,3)
    IGTAKA = COUNT (GTAKA(:))
!
    IF (IGTAKA > 0) THEN
      ZEW(:) = ZRAR(:) / 8.
      CALL ELEC_INI_NI_TAKAH(ZEW, ZDQRAR_SG, XTAKA_TM)
!
      WHERE (ZDQRAR_SG(:) < 0.)
        ZSAUNSK_SG(:) = XSKN_TAK
        ZSAUNSM_SG(:) = XSMN
        ZSAUNSN_SG(:) = XSNN
      ENDWHERE
    ENDIF
!
    DEALLOCATE( GTAKA )
  END IF
!
END SUBROUTINE ELEC_INI_NI_PROCESS
! 
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_INI_NI_SAP98(ZRAR, ZDQRAR_AUX)
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)    :: ZRAR
REAL, DIMENSION(:), INTENT(INOUT) :: ZDQRAR_AUX   ! q= f(RAR,T) in Saunders and
                                                  ! Peck's equation
!
  ZDQRAR_AUX(:) = 0.
!
! positive region : Mansell et al., 2005
  WHERE (GELEC(:,4) .AND. ZRAR(:) > ZRAR_CRIT(:))
    ZDQRAR_AUX(:) = MAX(0., 6.74 * (ZRAR(:) - ZRAR_CRIT(:)) * 1.E-15)
  ENDWHERE
!
! negative region : Mansell et al. 2005
  WHERE (GELEC(:,4) .AND. ZRAR(:) < ZRAR_CRIT(:))
    ZDQRAR_AUX(:) = MIN(0., 3.9 * (ZRAR_CRIT(:) - 0.1) *                   &
                           (4.0 * ((ZRAR(:) - (ZRAR_CRIT(:) + 0.1) / 2.) / &
                                  (ZRAR_CRIT(:) - 0.1))**2 - 1.) * 1.E-15)
  ENDWHERE
!
END SUBROUTINE ELEC_INI_NI_SAP98
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_INI_NI_SAUNQ(ZEW, ZDQLWC_AUX)
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN)    :: ZEW
REAL, DIMENSION(:), INTENT(INOUT) :: ZDQLWC_AUX  ! q= f(RAR or EW,T) in Saunders
                                                 !... equation
!
! For temperature lower than -30C and higher than -40C, value of q at -30C
!
  ALLOCATE ( IVEC1(IGSAUN) )
  ALLOCATE ( IVEC2(IGSAUN) )
  ALLOCATE ( ZVEC1(IGSAUN) )
  ALLOCATE ( ZVEC2(IGSAUN) )
  ALLOCATE ( ZDQLWC_OPT(IGSAUN) )
!
  ZDQLWC_OPT(:) = 0.
  IVEC1(:) = 0
  IVEC2(:) = 0
!
  ZVEC1(:) = PACK( ZZT(:), MASK=GSAUN(:))
  ZVEC2(:) = PACK( ZEW(:), MASK=GSAUN(:))
  ZDQLWC_OPT(:) = PACK( ZDQLWC_AUX(:), MASK=GSAUN )
!
! Temperature index (0C --> -40C)
  ZVEC1(1:IGSAUN) = MAX( 1.00001, MIN( REAL(NIND_TEMP)-0.00001, &
                                     (ZVEC1(1:IGSAUN) - XTT - 1.)/(-1.) ) )
  IVEC1(1:IGSAUN) = INT( ZVEC1(1:IGSAUN) )
  ZVEC1(1:IGSAUN) = ZVEC1(1:IGSAUN) - REAL(IVEC1(1:IGSAUN))
!
! LWC index (0.01 g.m^-3 --> 10 g.m^-3)
  WHERE (ZVEC2(:) >= 0.01 .AND. ZVEC2(:) < 0.1)
    ZVEC2(:) = MAX( 1.00001, MIN( REAL(10)-0.00001, &
                                  ZVEC2(:) * 100. ))
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
  WHERE (ZVEC2(:) >= 0.1 .AND. ZVEC2(:) < 1. .AND. IVEC2(:) == 0)
    ZVEC2(:) = MAX( 10.00001, MIN( REAL(19)-0.00001, &
                                   ZVEC2(:) * 10. + 9. ) )
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
  WHERE ((ZVEC2(:) >= 1.) .AND. ZVEC2(:) <= 10.)
    ZVEC2(:) = MAX( 19.00001, MIN( REAL(NIND_LWC)-0.00001, &
                                   ZVEC2(:) + 18. ) )
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
! Interpolate XSAUNDER
  ZDQLWC_OPT(:) = BI_LIN_INTP_V( XSAUNDER, IVEC2, IVEC1, ZVEC2, ZVEC1, &
                                  IGSAUN )
  ZDQLWC_AUX(:) = UNPACK( ZDQLWC_OPT(:), MASK=GSAUN, FIELD=0.0 )
!
  DEALLOCATE( IVEC1 )
  DEALLOCATE( IVEC2 )
  DEALLOCATE( ZVEC1 )
  DEALLOCATE( ZVEC2 )
  DEALLOCATE( ZDQLWC_OPT )
!
END SUBROUTINE ELEC_INI_NI_SAUNQ
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_INI_NI_TAKAH(ZEW, ZDQTAKA_AUX, XTAKA_AUX)
!
IMPLICIT NONE
!
REAL, DIMENSION(IMICRO) :: ZEW
REAL, DIMENSION(IMICRO) :: ZDQTAKA_AUX
REAL, DIMENSION(NIND_LWC+1,NIND_TEMP+1) :: XTAKA_AUX  !XMANSELL or XTAKA_TM)
!
!
  ALLOCATE ( IVEC1(IGTAKA) )
  ALLOCATE ( IVEC2(IGTAKA) )
  ALLOCATE ( ZVEC1(IGTAKA) )
  ALLOCATE ( ZVEC2(IGTAKA) )
  ALLOCATE ( ZDQTAKA_OPT(IGTAKA) )

  ZDQTAKA_OPT(:) = 0.
  IVEC1(:) = 0
  IVEC2(:) = 0
!
  ZVEC1(:) = PACK( ZZT(:), MASK=GTAKA )
  ZVEC2(:) = PACK( ZEW(:), MASK=GTAKA )
  ZDQTAKA_OPT(:) = PACK( ZDQTAKA_AUX(:), MASK=GTAKA )
!
! Temperature index (0C --> -40C)
  ZVEC1(1:IGTAKA) = MAX( 1.00001, MIN( REAL(NIND_TEMP)-0.00001, &
                                     (ZVEC1(1:IGTAKA) - XTT - 1.)/(-1.) ) )
  IVEC1(1:IGTAKA) = INT( ZVEC1(1:IGTAKA) )
  ZVEC1(1:IGTAKA) = ZVEC1(1:IGTAKA) - REAL(IVEC1(1:IGTAKA))
!
! LWC index (0.01 g.m^-3 --> 10 g.m^-3)
  WHERE (ZVEC2(:) >= 0.01 .AND. ZVEC2(:) < 0.1)
    ZVEC2(:) = MAX( 1.00001, MIN( REAL(10)-0.00001, &
                                  ZVEC2(:) * 100. ))
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
  WHERE (ZVEC2(:) >= 0.1 .AND. ZVEC2(:) < 1. .AND. IVEC2(:) == 0)
    ZVEC2(:) = MAX( 10.00001, MIN( REAL(19)-0.00001, &
                                   ZVEC2(:) * 10. + 9. ) )
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
  WHERE (ZVEC2(:) >= 1. .AND. ZVEC2(:) <= 10.)
    ZVEC2(:) = MAX( 19.00001, MIN( REAL(NIND_LWC)-0.00001, &
                                   ZVEC2(:) + 18. ) )
    IVEC2(:) = INT(ZVEC2(:))
    ZVEC2(:) = ZVEC2(:) - REAL(IVEC2(:))
  ENDWHERE
!
! Interpolate XMANSELL or XTAKA_TM
  ZDQTAKA_OPT(:) = BI_LIN_INTP_V( XTAKA_AUX, IVEC2, IVEC1, ZVEC2, ZVEC1, &
                                  IGTAKA )
  ZDQTAKA_AUX(:) = UNPACK( ZDQTAKA_OPT(:), MASK=GTAKA, FIELD=0.0 )
!
  DEALLOCATE( IVEC1 )
  DEALLOCATE( IVEC2 )
  DEALLOCATE( ZVEC1 )
  DEALLOCATE( ZVEC2 )
  DEALLOCATE( ZDQTAKA_OPT )
!
END SUBROUTINE ELEC_INI_NI_TAKAH
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_IAGGS_B()
!
!   Purpose : compute charge separation process during the collision
!             between ice and snow
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*      1. Collision efficiency
!
  ZCOLIS(:) = XCOLIS * EXP(XCOLEXIS * (ZZT(:) - XTT))
!
!*      2. Charging process following Helsdon and Farley (1987)
!
  IF (CNI_CHARGING == 'HELFA') THEN
    ZWQ1(:,7) = 0.
!
    WHERE (ZRIS(:) > XRTMIN_ELEC(4) .AND. ZCIT(:) > 0.0 .AND. &
           ZRIT(:) > XRTMIN_ELEC(4) .AND.                     &
           ZRST(:) > XRTMIN_ELEC(5))
      ZWQ1(:,7) = XFQIAGGSBH * ZZW(:) * ZCIT(:) / ZRIT(:)
      ZWQ1(:,7) = ZWQ1(:,7) * (1. - ZCOLIS(:)) / ZCOLIS(:) 
!
! Temperature dependance of the charge transferred
      ZWQ1(:,7) = ZWQ1(:,7) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
      ZWQ1(:,7) = ZWQ1(:,7) / ZRHODREF(:)
!
      ZQSS(:) = ZQSS(:) + ZWQ1(:,7)
      ZQIS(:) = ZQIS(:) - ZWQ1(:,7)
    END WHERE
  END IF
!
!*      3. Charging process following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    ZWQ1(:,7) = 0.
    WHERE (GELEC(:,1) .AND. ZDELTALWC(:) > 0. .AND.               &
           ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZRSS(:) > ZRSMIN_ELEC(5))
      ZWQ1(:,7) = XFQIAGGSBG  * (1 - ZCOLIS(:)) *          &
                  ZRHODREF(:)**(-4. * XCEXVT + 4. / XBI) * &
                  ZCIT(:)**(1 - 4. / XBI) *                &
                  ZDELTALWC(:) * ZFT(:) *                  &
                  ZLBDAS(:)**(XCXS - 2. - 4. * XDS) *      &
                  (XAI * MOMG(XALPHAI, XNUI, XBI) /        &
                  ZRIT(:))**(-4 / XBI)
!
! Dq is limited to XLIM_NI_IS
      ZLIMIT(:) = XLIM_NI_IS * ZZW(:) * ZCIT(:) * &
                  (1 - ZCOLIS(:)) / (ZRIT(:) * ZCOLIS(:))
      ZWQ1(:,7) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,7)) ), ZWQ1(:,7) )   
      ZWQ1(:,7) = ZWQ1(:,7) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,7) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,7) = ZWQ1(:,7) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
    ZQSS(:) = ZQSS(:) + ZWQ1(:,7)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,7)
  END IF
!
!*      4. Charging process based on EW: SAUN1/SAUN2, TEEWC
!*         following Saunders et al. (1991), Takahashi via Tsenova and Mitzeva (2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    ZWQ1(:,7) = 0.
!
    WHERE (GELEC(:,1) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZDQLWC(:) /= 0.)
      ZWQ1(:,7) = XFQIAGGSBS * (1 - ZCOLIS(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN(:)) *                          &
                  ZFQIAGGS(:) * ZDQLWC(:) *                               &
                  ZCIT(:)**(1 - ZSAUNIM(:) / XBI) *                       &
                  ZLBDAS(:)**(XCXS - 2.- XDS * (1. + ZSAUNIN(:))) *       &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM(:) / XBI)
!
! Dq is limited to XLIM_NI_IS
      ZLIMIT(:) = XLIM_NI_IS * ZZW(:) * ZCIT(:) * &
                  (1 - ZCOLIS(:)) / (ZRIT(:) * ZCOLIS(:))
      ZWQ1(:,7) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,7)) ), ZWQ1(:,7) )
      ZWQ1(:,7) = ZWQ1(:,7) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,7) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,7) = ZWQ1(:,7) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
    ZQSS(:) = ZQSS(:) + ZWQ1(:,7)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,7)   
!
  END IF
!
!*      5. Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*         following Saunders and Peck (1998) or
!*                   Brooks et al., 1997 (with/out anomalies) or
!*                   Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
!
    IF (CNI_CHARGING /= 'TERAR') THEN
       ZFQIAGGS(:)   = XFQIAGGSP
       WHERE (ZDQRAR_IS(:) < 0.)
         ZFQIAGGS(:) = XFQIAGGSN
       ENDWHERE
    ELSE
       ZFQIAGGS(:) = XFQIAGGSP_TAK
       WHERE (ZDQRAR_IS(:) <0.)
         ZFQIAGGS(:) = XFQIAGGSN_TAK
       ENDWHERE
    ENDIF
!
    ZWQ1(:,7) = 0.
!
    WHERE (GELEC(:,1) .AND. ZDQRAR_IS(:) /= 0. .AND. &
           ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZRSS(:) > ZRSMIN_ELEC(5))
      ZWQ1(:,7) = XFQIAGGSBS * (1 - ZCOLIS(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IS(:)) *                       &
                  ZFQIAGGS(:) * ZDQRAR_IS(:) *                            &
                  ZCIT(:)**(1 - ZSAUNIM_IS(:) / XBI) *                    &
                  ZLBDAS(:)**(XCXS - 2.- XDS * (1. + ZSAUNIN_IS(:))) *    &
                 (ZRHODREF(:) * ZRIT(:)/XAIGAMMABI)**(ZSAUNIM_IS(:) / XBI)
!
! Dq is limited to XLIM_NI_IS
      ZLIMIT(:) = XLIM_NI_IS * ZZW(:) * ZCIT(:) * &
                  (1 - ZCOLIS(:)) / (ZRIT(:) * ZCOLIS(:))
      ZWQ1(:,7) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,7)) ), ZWQ1(:,7) )
      ZWQ1(:,7) = ZWQ1(:,7) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,7) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,7) = ZWQ1(:,7) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
    ZQSS(:) = ZQSS(:) + ZWQ1(:,7)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,7)   
  END IF
!
!*      6. Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    ZWQ1(:,7) = 0.
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,1) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND.   &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZDQLWC(:) /= 0.)
      ZWQ1(:,7) = XFQIAGGSBT1 * (1.0 - ZCOLIS(:)) * ZRHOCOR(:) *          &
                  ZCIT(:) * ZLBDAS(:)**XCXS * ZDQLWC(:) *                 &
                  MIN( XFQIAGGSBT2 / (ZLBDAS(:)**(2. + XDS)) ,            &
                       XFQIAGGSBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./XBI) * &
                       ZRIT(:)**(2. / XBI) /                              &
                      (ZCIT(:)**(2. / XBI) * ZLBDAS(:)**(2. + 2. * XDS)))
!
! Dq is limited to XLIM_NI_IS
      ZLIMIT(:) = XLIM_NI_IS * ZZW(:) * ZCIT(:) * &
                  (1 - ZCOLIS(:)) / (ZRIT(:) * ZCOLIS(:))
      ZWQ1(:,7) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,7)) ), ZWQ1(:,7) )
      ZWQ1(:,7) = ZWQ1(:,7) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,7) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,7) = ZWQ1(:,7) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
    ZQSS(:) = ZQSS(:) + ZWQ1(:,7)
    ZQIS(:) = ZQIS(:) - ZWQ1(:,7)
  END IF
!
!
END SUBROUTINE ELEC_IAGGS_B
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ELEC_IDRYG_B()
!
!   Purpose : compute charge separation process during the dry collision
!             between ice and graupeln
!
!*	0.	DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE THE COLLECTION EFFICIENCY
!               ---------------------------------
!
  ZCOLIG(:) = XCOLIG * EXP(XCOLEXIG * (ZZT(:) - XTT))
!
!*      2.      COMPUTE THE CHARGE SEPARATION DURING IDRYG_BOUN
!               -----------------------------------------------
!
!*      2.1     Helsdon and  Farley (1987)
!
  IF (CNI_CHARGING == 'HELFA') THEN
    ZWQ1(:,3) = 0.
    WHERE (ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZCIT(:) > 0.0 .AND.            &
           ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6)) 
      ZWQ1(:,3) = XHIDRYG * ZZW1(:,2) * ZCIT(:) / ZRIT(:)    
      ZWQ1(:,3) = ZWQ1(:,3) * (1. - ZCOLIG(:)) / ZCOLIG(:)        ! QIDRYG_boun
!
! Temperature dependance of the charge transfered
      ZWQ1(:,3) = ZWQ1(:,3) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
      ZWQ1(:,3) = ZWQ1(:,3) / ZRHODREF(:)
    END WHERE
  END IF
!
!
!*      2.2     Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    ZWQ1(:,3) = 0.
!
    WHERE (GELEC(:,2) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND.        &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ZDELTALWC(:) > 0.) 
      ZWQ1(:,3) = XFQIDRYGBG * XLBQIDRYGBG * (1 - ZCOLIG) * &
                  ZRHODREF(:)**(-4. * XCEXVT + 4. / XBI) *  &
                  ZCIT(:)**(1 - 4. / XBI) *                 &
                  ZDELTALWC(:) * ZFT(:) *                   &
                  ZLBDAG(:)**(XCXG - 2. - 4. * XDG) *       &
                  (XAI * MOMG(XALPHAI, XNUI, XBI) /         &
                  ZRIT(:))**(-4 / XBI)
!
! Dq limited to XLIM_NI_IG
      ZLIMIT(:) = XLIM_NI_IG * ZZW1(:,2) * ZCIT(:) * (1 - ZCOLIG(:)) / &
                  (ZRIT(:) * ZCOLIG(:))
      ZWQ1(:,3) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,3)) ), ZWQ1(:,3) )
      ZWQ1(:,3) = ZWQ1(:,3) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
!
    WHERE (ZWQ1(:,3) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,3) = ZWQ1(:,3) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
  END IF
!
!
!*     2.3 Charging process based on EW: SAUN1/SAUN2, TEEWC
!*         following Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    ZWQ1(:,3) = 0.
!
    WHERE (GELEC(:,2) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ZDQLWC(:) /= 0.)
      ZWQ1(:,3) = XFQIDRYGBS * (1. - ZCOLIG(:)) *                    &
                  ZRHOCOR(:)**(1. + ZSAUNIN(:)) *                    &
                  ZFQIDRYGBS(:) * ZDQLWC(:) *                        &
                  ZCIT(:)**(1. - ZSAUNIM(:) / XBI) *                 &
                  ZLBDAG(:)**(XCXG - 2. - XDG * (1. + ZSAUNIN(:))) * &
                 (ZRHODREF(:) * ZRIT(:)/XAIGAMMABI)**(ZSAUNIM(:) / XBI)
!
! Dq is limited to XLIM_NI_IG
      ZLIMIT(:) = XLIM_NI_IG * ZZW1(:,2) * ZCIT(:) * (1 - ZCOLIG(:)) / &
                  (ZRIT(:) * ZCOLIG(:))
      ZWQ1(:,3) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,3)) ), ZWQ1(:,3) )
      ZWQ1(:,3) = ZWQ1(:,3) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,3) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,3) = ZWQ1(:,3) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
  END IF
!
!
!*     2.4 Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*         following Saunders and Peck (1998) or
!*                   Brooks et al., 1997 (with/out anomalies) or
!*                   Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
!
    IF (CNI_CHARGING /= 'TERAR') THEN
       ZFQIDRYGBS(:) = XFQIDRYGBSP
       WHERE (ZDQRAR_IG(:) < 0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN
       ENDWHERE
    ELSE
       ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
       WHERE (ZDQRAR_IG(:) <0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
       ENDWHERE
    END IF
!
    ZWQ1(:,3) = 0.
!
    WHERE (GELEC(:,2) .AND. ZDQRAR_IG(:) /= 0. .AND. &
           ZRIS(:) > ZRSMIN_ELEC(4) .AND. ZRGS(:) > ZRSMIN_ELEC(6)) 
      ZWQ1(:,3) = XFQIDRYGBS * (1. - ZCOLIG(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IG(:)) *                        &
                  ZFQIDRYGBS(:) * ZDQRAR_IG(:) *                           &
                  ZCIT(:)**(1 - ZSAUNIM_IG(:) / XBI) *                     &
                  ZLBDAG(:)**(XCXG - 2. - XDG * (1. + ZSAUNIN_IG(:))) *    &
                 (ZRHODREF(:) * ZRIT(:)/XAIGAMMABI)**(ZSAUNIM_IG(:) / XBI)
!
! Dq is limited to XLIM_NI_IG
      ZLIMIT(:) = XLIM_NI_IG * ZZW1(:,2) * ZCIT(:) * (1 - ZCOLIG(:)) / &
                  (ZRIT(:) * ZCOLIG(:))
      ZWQ1(:,3) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,3)) ), ZWQ1(:,3) )
      ZWQ1(:,3) = ZWQ1(:,3) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,3) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,3) = ZWQ1(:,3) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
  END IF
!
!
!*      2.5     Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    ZWQ1(:,3) = 0.
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,2) .AND. ZRIS(:) > ZRSMIN_ELEC(4) .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ZDQLWC(:) /= 0.)
      ZWQ1(:,3) = XFQIDRYGBT1 * (1. - ZCOLIG(:)) * ZRHOCOR(:) *           &
                  ZCIT(:) * ZLBDAG(:)**XCXG * ZDQLWC(:) *                 &
                  MIN( XFQIDRYGBT2 / (ZLBDAG(:)**(2. + XDG)),             &
                       XFQIDRYGBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./XBI) * &
                       ZRIT(:)**(2. / XBI) / (ZCIT(:)**(2. / XBI) *       &
                       ZLBDAG(:)**(2. + 2. * XDG)) )
!
! Dq is limited to XLIM_NI_IG
      ZLIMIT(:) = XLIM_NI_IG * ZZW1(:,2) * ZCIT(:) * (1 - ZCOLIG(:)) / &
                  (ZRIT(:) * ZCOLIG(:))
      ZWQ1(:,3) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,3)) ), ZWQ1(:,3) )
      ZWQ1(:,3) = ZWQ1(:,3) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,3) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,3) = ZWQ1(:,3) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
  END IF
!
!
END SUBROUTINE ELEC_IDRYG_B
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE  ELEC_SDRYG_B() 
!
!   Purpose : compute the charge separation during the dry collision
!             between snow and graupeln
!
!*	0.	DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE THE COLLECTION EFFICIENCY
!               ---------------------------------
!
  ZCOLSG(:) = XCOLSG * EXP (XCOLEXSG * (ZZT(:) - XTT))
!
!*      2.      COMPUTE THE CHARGE SEPARATION DURING SDRYG_BOUN
!               -----------------------------------------------
!
!*      2.1     Helsdon and Farley (1987)
! 
  IF (CNI_CHARGING == 'HELFA') THEN
    ZWQ1(:,5) = 0. 
!
    WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ZRSS(:) > ZRSMIN_ELEC(5) .AND. &
           ZLBDAS(:) > 0. .AND. ZLBDAG(:) > 0.)
      ZWQ1(:,5) = ZWQ1(:,10) * XFQSDRYGBH * ZRHODREF(:)**(-XCEXVT) *   &
                 (1. - ZCOLSG(:)) *                                    &
                  ZLBDAS(:)**(XCXS) * ZLBDAG(:)**(XCXG) *              &
                 (XLBQSDRYGB4H * ZLBDAS(:)**(-2.) +                    &
                  XLBQSDRYGB5H * ZLBDAS(:)**(-1.) * ZLBDAG(:)**(-1.) + &
                  XLBQSDRYGB6H * ZLBDAG(:)**(-2.))
!
! Temperature dependance of the charge transfered
      ZWQ1(:,5) = ZWQ1(:,5) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
      ZWQ1(:,5) = ZWQ1(:,5) / ZRHODREF(:)
    ENDWHERE
  ENDIF
!
!
!*      2.2     Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN 
    ZWQ1(:,5) = 0. 
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,3) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND.        &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZDELTALWC(:) > 0.) 
      ZWQ1(:,5) = XFQSDRYGBG * (1. - ZCOLSG(:)) *                      &
                  ZRHODREF(:)**(-4. * XCEXVT) *                        &
                  ZFT(:) * ZDELTALWC(:) *                              &
                  ZLBDAG(:)**XCXG * ZLBDAS(:)**XCXS *                  &
                 (XLBQSDRYGB4G * ZLBDAS(:)**(-4.) * ZLBDAG(:)**(-2.) + &
                  XLBQSDRYGB5G * ZLBDAS(:)**(-5.) * ZLBDAG(:)**(-1.) + &
                  XLBQSDRYGB6G * ZLBDAS(:)**(-6.)) *                   &
                  ZWQ1(:,10)
!
! Dq is limited to XLIM_NI_SG
      ZLIMIT(:) = XLIM_NI_SG * ZAUX1(:) * XAUX_LIM *                &
                  ZRHOCOR(:) * (1. - ZCOLSG(:)) *                   &
                  ZLBDAS(:)**(XCXS) * ZLBDAG(:)**(XCXG) *           &
                 (XAUX_LIM1 * ZLBDAS(:)**(-2.) +                    &
                  XAUX_LIM2 * ZLBDAS(:)**(-1.) * ZLBDAG(:)**(-1.) + &
                  XAUX_LIM3 * ZLBDAG(:)**(-2.))
      ZWQ1(:,5) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,5)) ), ZWQ1(:,5))
      ZWQ1(:,5) = ZWQ1(:,5) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,5) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,5) = ZWQ1(:,5) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
  END IF
!
!*     2.3 Charging process based on EW: SAUN1/SAUN2, TEEWC
!*         following Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
!
    ZWQ1(:,5) = 0.
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,3) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZDQLWC(:) /= 0.) 
!
!    ZWQ1(:,5) = ZWQ3(:) If graupel gains positive charge ZDQLWC(:) > 0.
!    ZWQ1(:,5) = ZWQ4(:) If graupel gains negative charge ZDQLWC(:) < 0.
      ZWQ1(:,5) = ZWQ3(:) * (0.5 + SIGN(0.5,ZDQLWC(:))) + &
                  ZWQ4(:) * (0.5 - SIGN(0.5,ZDQLWC(:)))
!
      ZWQ1(:,5) = ZWQ1(:,5) * XFQSDRYGBS * (1. - ZCOLSG(:)) *                 &
                  ZRHOCOR(:)**(1. + ZSAUNSN(:)) *                             &
                  ZSAUNSK(:) * ZDQLWC(:) *                                    &
                  ZLBDAG(:)**XCXG * ZLBDAS(:)**XCXS *                         &
                ( ZLBQSDRYGB1S(:) / (ZLBDAS(:)**ZSAUNSM(:) *ZLBDAG(:)**2)   + &
                  ZLBQSDRYGB2S(:) / (ZLBDAS(:)**( 1.+ZSAUNSM(:))*ZLBDAG(:)) + &
                  ZLBQSDRYGB3S(:) /  ZLBDAS(:)**(2.+ZSAUNSM(:)) )
!
! Dq is limited to XLIM_NI_SG
      ZLIMIT(:) = XLIM_NI_SG * ZAUX1(:) * XAUX_LIM *                &
                  ZRHOCOR(:) * (1. - ZCOLSG(:)) *                   &
                  ZLBDAS(:)**(XCXS) * ZLBDAG(:)**(XCXG) *           &
                ( XAUX_LIM1 / ZLBDAS(:)**2           +              &
                  XAUX_LIM2 /(ZLBDAS(:) * ZLBDAG(:)) +              &
                  XAUX_LIM3 / ZLBDAG(:)**2 )
      ZWQ1(:,5) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,5)) ), ZWQ1(:,5))
      ZWQ1(:,5) = ZWQ1(:,5) / ZRHODREF(:)
    ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,5) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,5) = ZWQ1(:,5) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
!
  END IF
!
!
!*     2.4 Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*         following Saunders and Peck (1998) or
!*                   Brooks et al., 1997 (with/out anomalies) or
!*                   Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
!
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
    WHERE (ZDQRAR_SG(:) < 0.)
      ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
      ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
      ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ENDWHERE
!
    ZWQ1(:,5) = 0.  
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,3) .AND. ZDQRAR_SG(:) /= 0. .AND. &
           ZRGS(:) > ZRSMIN_ELEC(6) .AND. ZRSS(:) > ZRSMIN_ELEC(5))
      ZWQ1(:,5) = ZWQ3(:) * (0.5+SIGN(0.5,ZDQRAR_SG(:))) + &
                  ZWQ4(:) * (0.5-SIGN(0.5,ZDQRAR_SG(:)))
!
      ZWQ1(:,5) = ZWQ1(:,5) * XFQSDRYGBS * (1. - ZCOLSG(:)) *                 &
                  ZRHOCOR(:)**(1. + ZSAUNSN_SG(:)) *                          &
                  ZSAUNSK_SG(:) * ZDQRAR_SG(:) *                              &
                  ZLBDAG(:)**XCXG * ZLBDAS(:)**XCXS *                         &
                 (ZLBQSDRYGB1S(:)/(ZLBDAS(:)**ZSAUNSM_SG(:) * ZLBDAG(:)**2) + &
                  ZLBQSDRYGB2S(:)/(ZLBDAS(:)**(1.+ZSAUNSM_SG(:))*ZLBDAG(:)) + &
                  ZLBQSDRYGB3S(:)/ ZLBDAS(:)**(2.+ZSAUNSM_SG(:)) )
!
!
! Dq is limited to XLIM_NI_SG
      ZLIMIT(:) = XLIM_NI_SG * ZAUX1(:) * XAUX_LIM *                &
                  ZRHOCOR(:) * (1. - ZCOLSG(:)) *                   &
                  ZLBDAS(:)**(XCXS) * ZLBDAG(:)**(XCXG) *           &
                ( XAUX_LIM1 / ZLBDAS(:)**2           +              &
                  XAUX_LIM2 /(ZLBDAS(:) * ZLBDAG(:)) +              &
                  XAUX_LIM3 / ZLBDAG(:)**2 )
      ZWQ1(:,5) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,5)) ), ZWQ1(:,5))
      ZWQ1(:,5) = ZWQ1(:,5) / ZRHODREF(:)
    ENDWHERE
!
! For temperature lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,5) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,5) = ZWQ1(:,5) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
  END IF
!
!
!*      2.5     Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN 
    ZWQ1(:,5) = 0.
    ZLIMIT(:) = 0.
!
    WHERE (GELEC(:,3) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND. &
           ZRSS(:) > ZRSMIN_ELEC(5) .AND. ZDQLWC(:) /= 0.)
      ZWQ1(:,5) = XFQSDRYGBT1 * (1. - ZCOLSG(:)) * ZRHOCOR(:) *             &
                  ZLBDAG(:)**XCXG * ZLBDAS(:)**XCXS * ZDQLWC(:) *           &
                  MIN(10. * (                                               &
                   ABS(XFQSDRYGBT2 / (ZLBDAG(:)**XDG * ZLBDAS(:)**2.) -     & 
                       XFQSDRYGBT3 / (ZLBDAS(:)**(2. + XDS))) +             &
                   ABS(XFQSDRYGBT4 / (ZLBDAG(:)**(2.+XDG)) -                &
                       XFQSDRYGBT5 / (ZLBDAS(:)**XDS * ZLBDAG(:)**2.)) +    &
                   ABS(XFQSDRYGBT6 / (ZLBDAG(:)**(1. + XDG) * ZLBDAS(:)) -  &
                       XFQSDRYGBT7 / (ZLBDAS(:)**(1. + XDS) * ZLBDAG(:)))), &
                   XFQSDRYGBT8 * ZRHOCOR(:) * ZWQ1(:,10) *                  &
                  (XFQSDRYGBT9  / (ZLBDAS(:)**2. * ZLBDAG(:)**2.) +         &
                   XFQSDRYGBT10 / (ZLBDAS(:)**4.) +                         &
                   XFQSDRYGBT11 / (ZLBDAS(:)**3. * ZLBDAG(:))))
!
! Dq is limited to XLIM_NI_SG
      ZLIMIT(:) = XLIM_NI_SG * ZAUX1(:) * XAUX_LIM *                &
                  ZRHOCOR(:) * (1. - ZCOLSG(:)) *                   &
                  ZLBDAS(:)**(XCXS) * ZLBDAG(:)**(XCXG) *           &
                ( XAUX_LIM1 / ZLBDAS(:)**2           +              &
                  XAUX_LIM2 /(ZLBDAS(:) * ZLBDAG(:)) +              &
                  XAUX_LIM3 / ZLBDAG(:)**2 )
      ZWQ1(:,5) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ1(:,5)) ), ZWQ1(:,5))
      ZWQ1(:,5) = ZWQ1(:,5) / ZRHODREF(:)
    ENDWHERE
!
! For temperature lower than -30C --> linear interpolation
    WHERE (ZWQ1(:,5) /= 0. .AND. ZZT(:) < (XTT-30.) .AND. ZZT(:) >= (XTT-40.))
      ZWQ1(:,5) = ZWQ1(:,5) * (ZZT(:) - XTT + 40.) / 10.
    ENDWHERE
  END IF
!
!
END SUBROUTINE ELEC_SDRYG_B
!
!------------------------------------------------------------------------------
!
  SUBROUTINE INDUCTIVE_PROCESS
!
! Computation of the charge transfer rate during inductive mechanism
! Only the bouncing droplet-graupel collision when the graupel is in the dry 
! growth mode is considered
! The electric field is limited to 100 kV/m
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!*      1. COMPUTE THE CHARGING RATE
!          -------------------------
!
  ZRATE_IND(:) = 0.
!
  WHERE (GIND(:) .AND.                                                 &
         ZEFIELDW(:) /= 0. .AND. ABS(ZEGS(:)) > XEGMIN .AND.           &
         ZLBDAG(:) > 0. .AND.                                          &
         ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRGS(:) > ZRSMIN_ELEC(6) .AND. &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRCS(:) > ZRSMIN_ELEC(2))
    ZRATE_IND(:) = XIND1 * ZLBDAG(:)**XCXG * ZRHOCOR(:) *                     &
                  (XIND2 * SIGN(MIN(100.E3, ABS(ZEFIELDW(:))), ZEFIELDW(:)) * &
                   ZLBDAG(:) **(-2.-XDG) -                                    &
                   XIND3 * ZEGS(:) * ZLBDAG(:)**(-XFG-XDG))
    ZRATE_IND(:) = ZRATE_IND(:) / ZRHODREF(:)
    ZQGS(:) = ZQGS(:) + ZRATE_IND(:)
    ZQCS(:) = ZQCS(:) - ZRATE_IND(:)
  END WHERE
!
END SUBROUTINE INDUCTIVE_PROCESS
!
!------------------------------------------------------------------------------
!
!
  FUNCTION BI_LIN_INTP_V(ZT, KI, KJ, PDX, PDY, KN)  RESULT(Y)
!
!                 |                   |
! ZT(KI(1),KJ(2))-|-------------------|-ZT(KI(2),KJ(2))
!                 |                   | 
!                 |                   |
!              x2-|-------|y(x1,x2)   |
!                 |       |           |
!              PDY|       |           |
!                 |       |           |
!                 |       |           |
!ZT( KI(1),KJ(1))-|-------------------|-ZT(KI(2),KJ(1))
!                 |  PDX  |x1         |
!                 |                   |
!
!*	0.	DECLARATIONS
!          	------------
!
IMPLICIT NONE
!
!*	0.2	Declaration of local variables
!
INTEGER                          :: KN        ! Size of the result vector
INTEGER,          DIMENSION(KN)  :: KI        ! Tabulated  coordinate
INTEGER,          DIMENSION(KN)  :: KJ        ! Tabulated  coordinate
REAL, INTENT(IN), DIMENSION(:,:) :: ZT        ! Tabulated data
REAL, INTENT(IN), DIMENSION(KN)  :: PDX, PDY  ! 
REAL,             DIMENSION(KN)  :: Y         ! Interpolated value
!
INTEGER                          :: JJ        ! Loop index
!
!*	1.	INTERPOLATION
!		-------------
!  
DO JJ = 1, KN
  Y(JJ) = (1.0 - PDX(JJ)) * (1.0 - PDY(JJ)) * ZT(KI(JJ),  KJ(JJ))   + &
           PDX(JJ)        * (1.0 - PDY(JJ)) * ZT(KI(JJ)+1,KJ(JJ))   + &
           PDX(JJ)        * PDY(JJ)         * ZT(KI(JJ)+1,KJ(JJ)+1) + &
          (1.0 - PDX(JJ)) * PDY(JJ)         * ZT(KI(JJ)  ,KJ(JJ)+1)
ENDDO
!
END FUNCTION BI_LIN_INTP_V
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RAIN_ICE_ELEC
