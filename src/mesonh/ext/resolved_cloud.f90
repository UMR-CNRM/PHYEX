!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_RESOLVED_CLOUD
!     ##########################
INTERFACE
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HELEC, HACTCCN, HSCONV, HMF_CLOUD,           &
                                  KRR, KSPLITR, KSPLITG, KMI, KTCOUNT,                 &
                                  HLBCX, HLBCY, TPFILE, HRAD, HTURBDIM,                &
                                  OSUBG_COND, OSIGMAS, HSUBG_AUCV,                     &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,              &
                                  PPABST, PTHT, PRT, PSIGS, PSIGQSAT, PMFCONV,         &
                                  PTHM, PRCM, PPABSTT,                                 &
                                  PW_ACT,PDTHRAD, PTHS, PRS, PSVT, PSVS, PSRCS, PCLDFR,&
                                  PICEFR,                                              &
                                  PCIT, OSEDIC, OACTIT, OSEDC, OSEDI,                  &
                                  ORAIN, OWARM, OHHONI, OCONVHG,                       &
                                  PCF_MF,PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,             &
                                  PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF,  &
                                  PINPRC,PINPRC3D,PINPRR,PINPRR3D, PEVAP3D,            &
                                  PINPRS,PINPRS3D,PINPRG,PINPRG3D,PINPRH,PINPRH3D,     &
                                  PSOLORG,PMI,                                         &
                                  PSPEEDC, PSPEEDR, PSPEEDS, PSPEEDG, PSPEEDH,         &
                                  PINDEP, PSUPSAT,  PNACT, PNPRO,PSSPRO, PRAINFR,      &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,              &
                                  PSEA,PTOWN          )   
!
USE MODD_IO, ONLY: TFILEDATA
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud
CHARACTER(LEN=4),         INTENT(IN)   :: HELEC    ! kind of electrical scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HACTCCN  ! kind of CCN activation scheme
                                                   ! paramerization
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KSPLITG  ! Number of small time step
                                       ! integrations for  ice  sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
TYPE(TFILEDATA),          INTENT(INOUT)   :: TPFILE   ! Output file
CHARACTER(len=4),         INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER(len=4),         INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL,                     INTENT(IN)   :: PSIGQSAT! coeff applied to qsat variance contribution
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABSTT  ! Pressure time t+Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRCM    ! Cloud water m.r. at time t-Dt
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_ACT ! W for CCN activation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS   ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS    ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT   ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS   ! Scalar variable sources
!
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
LOGICAL,                  INTENT(IN)    :: OSEDIC! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for ICE3            
LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                 ! activation through temp.
                                                 ! evolution in C2R2 and KHKO
LOGICAL,                  INTENT(IN)    :: OSEDC ! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for C2R2 or KHKO
LOGICAL,                  INTENT(IN)    :: OSEDI ! Switch to activate the
                                                 ! cloud crystal sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN ! Switch to activate the
                                                 ! raindrop formation
LOGICAL,                  INTENT(IN)    :: OWARM ! Control of the rain formation
                                                 !  by slow warm microphysical
                                                 !         processes
LOGICAL,                  INTENT(IN)    :: OHHONI! enable haze freezing
LOGICAL,                  INTENT(IN)    :: OCONVHG! Switch for conversion from
                                                  ! hail to graupel
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PWEIGHT_MF_CLOUD ! weight coefficient for the mass-flux cloud
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC   ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH   ! Hail instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRC3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRS3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRG3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRH3D ! sed flux of precip
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG  ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDC ! Cloud sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDR ! Rain sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDS ! Snow sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDG ! Graupel sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDH ! Hail sedimentation speed
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP   ! Cloud instant deposition
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT    !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO    !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   !sursat
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR  ! Rain fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLC_HRC !HighLow liquid content
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLC_HCF !HighLow liquid cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLI_HRI !HighLow ice content
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLI_HCF !HighLow ice clous fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA      ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN      ! Town fraction
!
END SUBROUTINE RESOLVED_CLOUD
END INTERFACE
END MODULE MODI_RESOLVED_CLOUD
!
!     ##################################################################################
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HELEC, HACTCCN, HSCONV, HMF_CLOUD,           &
                                  KRR, KSPLITR, KSPLITG, KMI, KTCOUNT,                 &
                                  HLBCX, HLBCY, TPFILE, HRAD, HTURBDIM,                &
                                  OSUBG_COND, OSIGMAS, HSUBG_AUCV,                     &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,              &
                                  PPABST, PTHT, PRT, PSIGS, PSIGQSAT, PMFCONV,         &
                                  PTHM, PRCM, PPABSTT,                                 &
                                  PW_ACT,PDTHRAD, PTHS, PRS, PSVT, PSVS, PSRCS, PCLDFR,&
                                  PICEFR,                                              &
                                  PCIT, OSEDIC, OACTIT, OSEDC, OSEDI,                  &
                                  ORAIN, OWARM, OHHONI, OCONVHG,                       &
                                  PCF_MF,PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,             &
                                  PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF,  &
                                  PINPRC,PINPRC3D,PINPRR,PINPRR3D, PEVAP3D,            &
                                  PINPRS,PINPRS3D,PINPRG,PINPRG3D,PINPRH,PINPRH3D,     &
                                  PSOLORG,PMI,                                         &
                                  PSPEEDC, PSPEEDR, PSPEEDS, PSPEEDG, PSPEEDH,         &
                                  PINDEP, PSUPSAT,  PNACT, PNPRO,PSSPRO, PRAINFR,      &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,              &
                                  PSEA,PTOWN          )   
!     ##################################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation
!!
!!
!!**  METHOD
!!    ------
!!      The main actions of this routine is to call the routines computing the
!!    microphysical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical routines, the physical tendencies are
!!    switched back to prognostic variables.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine SLOW_TERMS: Computes the explicit microphysical sources
!!      Subroutine FAST_TERMS: Performs the saturation adjustment for l
!!      Subroutine RAIN_ICE  : Computes the explicit microphysical sources for i
!!      Subroutine ICE_ADJUST: Performs the saturation adjustment for i+l
!!      MIN_ll,SUM3D_ll : distributed functions equivalent to MIN and SUM
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          CST%XP00               ! Reference pressure
!!          CST%XRD                ! Gaz  constant for dry air
!!          CST%XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and book2 of documentation ( routine RESOLVED_CLOUD )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/12/94
!!      Modifications: June 8, 1995 ( J.Stein )
!!                                   Cleaning to improve efficienty and clarity
!!                                  in agreement with the MESO-NH coding norm
!!                     March 1, 1996 ( J.Stein )
!!                                   store the cloud fraction
!!                     March 18, 1996 ( J.Stein )
!!                                   check that ZMASSPOS /= 0
!!                     Oct.  12, 1996 ( J.Stein )
!!                                   remove the negative values correction
!!                                   for the KES2 case
!!      Modifications: Dec 14, 1995 (J.-P. Pinty)
!!                                   Add the mixed-phase option
!!      Modifications: Jul 01, 1996 (J.-P. Pinty)
!!                                   Change arg. list in routine FAST_TERMS
!!      Modifications: Jan 27, 1997 (J.-P. Pinty)
!!                                   add W and SV in arg. list
!!      Modifications: March 23, 98 (E.Richard)
!!                                   correction of negative value based on
!!                                  rv+rc+ri and thetal or thetail conservation
!!      Modifications: April 08, 98 (J.-P. Lafore and V. Ducrocq )
!!                                  modify the  correction of negative values
!!      Modifications: June 08, 00  (J.-P. Pinty and J.-M. Cohard)
!!                                  add the C2R2 scheme
!!      Modifications: April 08, 01  (J.-P. Pinty)
!!                                  add the C3R5 scheme
!!      Modifications: July  21, 01  (J.-P. Pinty)
!!                                  Add OHHONI and PW_ACT (for haze freezing)
!!      Modifications: Sept 21, 01  (J.-P. Pinty)
!!                                  Add XCONC_CCN limitation
!!      Modifications: Nov  21, 02  (J.-P. Pinty)
!!                                  Add ICE4 and C3R5 options
!!                     June, 2005   (V. Masson)
!!                                  Technical change in interface for scalar arguments
!!      Modifications : March, 2006 (O.Geoffroy)
!!                                  Add KHKO scheme
!!      Modifications : March 2013  (O.Thouron)
!!                                  Add prognostic supersaturation
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      M.Mazoyer : 04/2016 : Temperature radiative tendency used for  
!!                            activation by cooling (OACTIT)
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!                     10/2016 M.Mazoyer New KHKO output fields
!!                    10/2016 (C.Lac) Add droplet deposition
!!      S.Riette  : 11/2016 : ice_adjust before and after rain_ice
!!                            ICE3/ICE4 modified, old version under LRED=F   
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 01/02/2019: ZRSMIN is now allocatable (instead of size of XRTMIN which was sometimes not allocated)
!  C. Lac         02/2019: add rain fraction as an output field
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  B. Vie         03/2020: LIMA negativity checks after turbulence, advection and microphysics budgets
!  B. Vie      03/03/2020: use DTHRAD instead of dT/dt in Smax diagnostic computation
!  P. Wautelet 11/06/2020: bugfix: correct ZSVS array indices
!  P. Wautelet 11/06/2020: bugfix: add "Non local correction for precipitating species" for ICE4
!  P. Wautelet + Benoit Vié 06/2020: improve removal of negative scalar variables + adapt the corresponding budgets
!  P. Wautelet 23/06/2020: remove ZSVS and ZSVT to improve code readability
!  P. Wautelet 30/06/2020: move removal of negative scalar variables to Sources_neg_correct
!  P. Wautelet 30/06/2020: remove non-local corrections
!  B. Vie         06/2020: add prognostic supersaturation for LIMA
!  C. Barthe   20/03/2023: to avoid duplicating sources, cloud electrification is integrated in the microphysics
!                          CELLS can be used with rain_ice with LRED=T and with LIMA with LPTSPLIT=T
!                          the adjustement for cloud electricity is also externalized
!  A. Marcel Jan 2025: bi-Gaussian PDF and associated subgrid precipitation
!  A. Marcel Jan 2025: relaxation of the small fraction assumption
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_BUDGET,           ONLY: TBUDGETS, TBUDGETS_PTR, TBUCONF
USE MODD_CH_AEROSOL,       ONLY: LORILAM, NSP, NCARB, NSOA
USE MODD_CST,              ONLY: CST
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_DUST,             ONLY: LDUST
USE MODD_ELEC_DESCR,       ONLY: ELEC_DESCR, LSEDIM_BEARD, LIAGGS_LATHAM
USE MODD_ELEC_n,           ONLY: XEFIELDU, XEFIELDV, XEFIELDW
USE MODD_ELEC_PARAM,       ONLY: ELEC_PARAM
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_NEB_n,            ONLY: NEBN, CCONDENS, CLAMBDA3
USE MODD_NSV,              ONLY: NSV, NSV_C1R3END, NSV_C2R2BEG, NSV_C2R2END,                       &
                                 NSV_LIMA_BEG, NSV_LIMA_END, NSV_LIMA_CCN_FREE, NSV_LIMA_IFN_FREE, &
                                 NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_NR,                            &
                                 NSV_AEREND, NSV_DSTEND, NSV_SLTEND,                               &
                                 NSV_ELECBEG, NSV_ELECEND, TNSV
USE MODD_PARAM_C2R2,       ONLY: LSUPSAT
USE MODD_PARAMETERS,       ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_ICE_n,      ONLY: CSEDIM, LADJ_BEFORE, LADJ_AFTER, LRED, PARAM_ICEN
USE MODD_PARAM_LIMA,       ONLY: LADJ, LPTSPLIT, LSPRO, NMOD_CCN, NMOD_IFN, NMOD_IMM, NMOM_I, PARAM_LIMA
USE MODD_PARAM_LIMA_WARM,  ONLY: PARAM_LIMA_WARM
USE MODD_PARAM_LIMA_COLD,  ONLY: PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED, ONLY: PARAM_LIMA_MIXED
USE MODD_RAIN_ICE_DESCR_n, ONLY: XRTMIN, RAIN_ICE_DESCRN
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAMN
USE MODD_REF,              ONLY: XTHVREFZ
USE MODD_SALT,             ONLY: LSALT
USE MODD_TURB_n,           ONLY: TURBN
!
USE MODE_ll
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
USE MODE_MPPDB
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,   ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_sources_neg_correct, only: Sources_neg_correct
!
USE MODI_AER2LIMA
USE MODI_C2R2_ADJUST
USE MODI_ELEC_ADJUST
USE MODI_FAST_TERMS
USE MODI_GET_HALO
USE MODI_ICE_ADJUST
USE MODI_ICE_ADJUST_ELEC
USE MODI_ION_SOURCE_ELEC
USE MODI_KHKO_NOTADJUST
USE MODI_LIMA
USE MODI_LIMA_ADJUST
USE MODI_LIMA_ADJUST_SPLIT
USE MODI_LIMA_COLD
USE MODI_LIMA_MIXED
USE MODI_LIMA_NOTADJUST
USE MODI_LIMA_WARM
USE MODI_RAIN_C2R2_KHKO
USE MODI_RAIN_ICE
USE MODI_RAIN_ICE_ELEC
USE MODI_RAIN_ICE_OLD
USE MODI_SHUMAN
#ifdef MNH_OPENACC
USE MODI_SHUMAN_DEVICE
#endif
USE MODI_SLOW_TERMS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud paramerization
CHARACTER(LEN=4),         INTENT(IN)   :: HELEC    ! kind of electrical scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HACTCCN  ! kind of CCN activation scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KSPLITG  ! Number of small time step
                                       ! integrations for  ice  sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
TYPE(TFILEDATA),          INTENT(INOUT)   :: TPFILE   ! Output file
CHARACTER(len=4),         INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER(len=4),         INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL,                     INTENT(IN)   :: PSIGQSAT! coeff applied to qsat variance contribution
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABSTT ! Pressure time t+Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRCM    ! Cloud water m.r. at time t-Dt
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_ACT ! W for CCN activation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS   ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS    ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT   ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS   ! Scalar variable sources
!
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
LOGICAL,                  INTENT(IN)    :: OSEDIC! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for ICE3            
LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                 ! activation through temp.
                                                 ! evolution in C2R2 and KHKO
LOGICAL,                  INTENT(IN)    :: OSEDC ! Switch to activate the
                                                 ! cloud droplet sedimentation
LOGICAL,                  INTENT(IN)    :: OSEDI ! Switch to activate the
                                                 ! cloud crystal sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN ! Switch to activate the
                                                 ! raindrop formation
LOGICAL,                  INTENT(IN)    :: OWARM ! Control of the rain formation
                                                 !  by slow warm microphysical
                                                 !         processes
LOGICAL,                  INTENT(IN)    :: OHHONI! enable haze freezing
LOGICAL,                  INTENT(IN)    :: OCONVHG! Switch for conversion from
                                                  ! hail to graupel
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PWEIGHT_MF_CLOUD
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC   ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH   ! Hail instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRC3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRS3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRG3D ! sed flux of precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRH3D ! sed flux of precip
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG  ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDC ! Cloud sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDR ! Rain sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDS ! Snow sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDG ! Graupel sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDH ! Hail sedimentation speed
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP   ! Cloud instant deposition
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT    !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO    !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   !sursat
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR  ! Rain fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLC_HRC !HighLow liquid content
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLC_HCF !HighLow liquid cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLI_HRI !HighLow ice content
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PHLI_HCF !HighLow ice clous fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA      ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN      ! Town fraction
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR,JSV       ! Loop index for the moist and scalar variables
INTEGER :: IIB           !  Define the physical domain
INTEGER :: IIE           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: IKU
INTEGER :: IINFO_ll      ! return code of parallel routine
INTEGER :: JK,JI,JL
!
!
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZDZZ
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZEXN
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZZZ
                                    ! model layer height
! REAL  :: ZMASSTOT                   ! total mass  for one water category
!                                     ! including the negative values
! REAL  :: ZMASSPOS                   ! total mass  for one water category
!                                     ! after removing the negative values
! REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
!
INTEGER                               :: ISVBEG ! first scalar index for microphysics
INTEGER                               :: ISVEND ! last  scalar index for microphysics
!UPG*PT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSVT   ! scalar variable for microphysics only
!UPG*PT

REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3), KRR) :: ZFPR
!
INTEGER                               :: JMOD, JMOD_IFN
LOGICAL                               :: GWEST,GEAST,GNORTH,GSOUTH
LOGICAL                               :: LMFCONV ! =SIZE(PMFCONV)!=0
! BVIE work array waiting for PINPRI
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)):: ZINPRI
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSIGQSAT2D
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZDUM
!
! variables for cloud electricity
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCND, ZDEP
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRCS_BEF, ZRIS_BEF
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQHT, ZQPIT, ZQNIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQHS, ZQPIS, ZQNIS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLATHAM_IAGGS ! E Function to simulate
                                                     ! enhancement of IAGGS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEFIELDW
LOGICAL :: GELEC ! if true, cloud electrification is activated
LOGICAL :: LLDTHRAD ! true if DTHRAD is used by LIMA
INTEGER  :: JIU,JJU,JKU
INTEGER  :: JII,JJI,JKI
!
ZSIGQSAT2D(:,:) = PSIGQSAT
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
IKU=SIZE(PZZ,3)
!
JIU =  size(PZZ, 1 )
JJU =  size(PZZ, 2 )
JKU =  size(PZZ, 3 )
!
CALL FILL_DIMPHYEX(YLDIMPHYEX, SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3))
!
GWEST  = ( HLBCX(1) /= 'CYCL' .AND. LWEST_ll() )
GEAST  = ( HLBCX(2) /= 'CYCL' .AND. LEAST_ll() )
GSOUTH = ( HLBCY(1) /= 'CYCL' .AND. LSOUTH_ll() )
GNORTH = ( HLBCY(2) /= 'CYCL' .AND. LNORTH_ll() )
!
LMFCONV=(SIZE(PMFCONV)/=0)
!
IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') THEN
  ISVBEG = NSV_C2R2BEG
  ISVEND = NSV_C2R2END
ELSE IF (HCLOUD == 'C3R5') THEN
  ISVBEG = NSV_C2R2BEG
  ISVEND = NSV_C1R3END
ELSE IF (HCLOUD == 'LIMA') THEN
  ISVBEG = NSV_LIMA_BEG
  IF (.NOT. LDUST .AND. .NOT. LSALT .AND. .NOT. LORILAM) THEN
    ISVEND = NSV_LIMA_END
  ELSE
    IF (LORILAM) THEN
      ISVEND = NSV_AEREND
    END IF
    IF (LDUST) THEN
      ISVEND = NSV_DSTEND
    END IF
    IF (LSALT) THEN
      ISVEND = NSV_SLTEND
    END IF
  END IF
ELSE
  ISVBEG = 0
  ISVEND = 0
END IF
!
!
!
!*        1.    From ORILAM to LIMA: 
!
IF (HCLOUD == 'LIMA' .AND. ((LORILAM).OR.(LDUST).OR.(LSALT))) THEN
! ORILAM : tendance s --> variable instant t
  ALLOCATE(ZSVT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NSV))
  DO JSV = 1, NSV
    ZSVT(:,:,:,JSV) = PSVS(:,:,:,JSV) * PTSTEP / PRHODJ(:,:,:)
  END DO

  CALL AER2LIMA(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), &
                PRT(IIB:IIE,IJB:IJE,IKB:IKE,1),&
                PPABST(IIB:IIE,IJB:IJE,IKB:IKE),&
                PTHT(IIB:IIE,IJB:IJE,IKB:IKE), &
                PZZ(IIB:IIE,IJB:IJE,IKB:IKE))

! LIMA : variable instant t --> tendance s
  PSVS(:,:,:,NSV_LIMA_CCN_FREE)   = ZSVT(:,:,:,NSV_LIMA_CCN_FREE) * &
                                    PRHODJ(:,:,:) / PTSTEP
  PSVS(:,:,:,NSV_LIMA_CCN_FREE+1) = ZSVT(:,:,:,NSV_LIMA_CCN_FREE+1) * &
                                    PRHODJ(:,:,:) / PTSTEP
  PSVS(:,:,:,NSV_LIMA_CCN_FREE+2) = ZSVT(:,:,:,NSV_LIMA_CCN_FREE+2) * &
                                    PRHODJ(:,:,:) / PTSTEP

  PSVS(:,:,:,NSV_LIMA_IFN_FREE)   = ZSVT(:,:,:,NSV_LIMA_IFN_FREE) * &
                                    PRHODJ(:,:,:) / PTSTEP
  PSVS(:,:,:,NSV_LIMA_IFN_FREE+1) = ZSVT(:,:,:,NSV_LIMA_IFN_FREE+1) * &
                                    PRHODJ(:,:,:) / PTSTEP
  !
  DEALLOCATE(ZSVT)
END IF

!UPG*PT
!
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!$acc kernels
PTHS(:,:,:) = PTHS(:,:,:) / PRHODJ(:,:,:)
!$acc loop independent
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO
!$acc end kernels
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
!$acc kernels
!$acc loop independent
  DO JSV = ISVBEG, ISVEND
    PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) / PRHODJ(:,:,:)
  ENDDO
!$acc end kernels
ENDIF
!
!  complete the lateral boundaries to avoid possible problems
!
!$acc kernels
!$acc loop independent
DO JI=1,JPHEXT
  PTHS(JI,:,:) = PTHS(IIB,:,:)
  PTHS(IIE+JI,:,:) = PTHS(IIE,:,:)
  PTHS(:,JI,:) = PTHS(:,IJB,:)
  PTHS(:,IJE+JI,:) = PTHS(:,IJE,:)
!
  PRS(JI,:,:,:) = PRS(IIB,:,:,:)
  PRS(IIE+JI,:,:,:) = PRS(IIE,:,:,:)
  PRS(:,JI,:,:) = PRS(:,IJB,:,:)
  PRS(:,IJE+JI,:,:) = PRS(:,IJE,:,:)
END DO
!
!  complete the physical boundaries to avoid some computations
!
IF( GWEST  )  PRT(:IIB-1,:,:,2:) = 0.0
IF( GEAST  )  PRT(IIE+1:,:,:,2:) = 0.0
IF( GSOUTH )  PRT(:,:IJB-1,:,2:) = 0.0
IF( GNORTH )  PRT(:,IJE+1:,:,2:) = 0.0
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
!$acc loop independent
DO JI=1,JPHEXT
  PSVS(JI,     :,      :, ISVBEG:ISVEND) = PSVS(IIB, :,   :, ISVBEG:ISVEND)
  PSVS(IIE+JI, :,      :, ISVBEG:ISVEND) = PSVS(IIE, :,   :, ISVBEG:ISVEND)
  PSVS(:,      JI,     :, ISVBEG:ISVEND) = PSVS(:,   IJB, :, ISVBEG:ISVEND)
  PSVS(:,      IJE+JI, :, ISVBEG:ISVEND) = PSVS(:,   IJE, :, ISVBEG:ISVEND)
END DO
 !
!  complete the physical boundaries to avoid some computations
!
  IF( GWEST  ) PSVT(:IIB-1, :,      :, ISVBEG:ISVEND) = 0.0
  IF( GEAST  ) PSVT(IIE+1:, :,      :, ISVBEG:ISVEND) = 0.0
  IF( GSOUTH ) PSVT(:,      :IJB-1, :, ISVBEG:ISVEND) = 0.0
  IF( GNORTH ) PSVT(:,      IJE+1:, :, ISVBEG:ISVEND) = 0.0
ENDIF
!
!  complete the vertical boundaries
!
PTHS(:,:,IKB-1) = PTHS(:,:,IKB)
PTHS(:,:,IKE+1) = PTHS(:,:,IKE)
!
PRS(:,:,IKB-1,:) = PRS(:,:,IKB,:)
PRS(:,:,IKE+1,:) = PRS(:,:,IKE,:)
!
PRT(:,:,IKB-1,:) = PRT(:,:,IKB,:)
PRT(:,:,IKE+1,:) = PRT(:,:,IKE,:)
!
IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO' &
                                           .OR. HCLOUD == 'LIMA') THEN
  PSVS(:,:,IKB-1,ISVBEG:ISVEND) = PSVS(:,:,IKB,ISVBEG:ISVEND)
  PSVS(:,:,IKE+1,ISVBEG:ISVEND) = PSVS(:,:,IKE,ISVBEG:ISVEND)
  PSVT(:,:,IKB-1,ISVBEG:ISVEND) = PSVT(:,:,IKB,ISVBEG:ISVEND)
  PSVT(:,:,IKE+1,ISVBEG:ISVEND) = PSVT(:,:,IKE,ISVBEG:ISVEND)
ENDIF
!$acc end kernels
!
! Same thing for cloud electricity
IF (HELEC(1:3) == 'ELE') THEN
  ! Transformation into physical tendencies
  DO JSV = NSV_ELECBEG, NSV_ELECEND
    PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) / PRHODJ(:,:,:)
  ENDDO
  !
  ! complete the lateral boundaries to avoid possible problems
  DO JI = 1, JPHEXT
    ! positive ion source
    PSVS(JI,:,:,NSV_ELECBEG)     = PSVS(IIB,:,:,NSV_ELECBEG)
    PSVS(IIE+JI,:,:,NSV_ELECBEG) = PSVS(IIE,:,:,NSV_ELECBEG)
    PSVS(:,JI,:,NSV_ELECBEG)     = PSVS(:,IJB,:,NSV_ELECBEG)
    PSVS(:,IJE+JI,:,NSV_ELECBEG) = PSVS(:,IJE,:,NSV_ELECBEG)
    ! source of hydrometeor charge
    PSVS(JI,:,:,NSV_ELECBEG+1:NSV_ELECEND-1)     = 0.0
    PSVS(IIE+JI,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
    PSVS(:,JI,:,NSV_ELECBEG+1:NSV_ELECEND-1)     = 0.0
    PSVS(:,IJE+JI,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
    ! negative ion source
    PSVS(JI,:,:,NSV_ELECEND)     = PSVS(IIB,:,:,NSV_ELECEND)
    PSVS(IIE+JI,:,:,NSV_ELECEND) = PSVS(IIE,:,:,NSV_ELECEND)
    PSVS(:,JI,:,NSV_ELECEND)     = PSVS(:,IJB,:,NSV_ELECEND)
    PSVS(:,IJE+JI,:,NSV_ELECEND) = PSVS(:,IJE,:,NSV_ELECEND)
  END DO
  !
  ! complete the physical boundaries to avoid some computations
  IF( GWEST  ) PSVT(IIB-1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  IF( GEAST  ) PSVT(IIE+1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  IF( GSOUTH ) PSVT(:,IJB-1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  IF( GNORTH ) PSVT(:,IJE+1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  !
  ! complete the vertical boundaries
  PSVS(:,:,IKB-1,NSV_ELECBEG) = PSVS(:,:,IKB,NSV_ELECBEG)    ! Positive ion
  PSVT(:,:,IKB-1,NSV_ELECBEG) = PSVT(:,:,IKB,NSV_ELECBEG)
  PSVS(:,:,IKB-1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0          ! Hydrometeor charge
  PSVS(:,:,IKE+1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  PSVT(:,:,IKB-1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  PSVT(:,:,IKE+1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
  PSVS(:,:,IKB-1,NSV_ELECEND) = PSVS(:,:,IKB,NSV_ELECEND)    ! Negative ion
  PSVT(:,:,IKB-1,NSV_ELECEND) = PSVT(:,:,IKB,NSV_ELECEND)
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
call Sources_neg_correct( hcloud, helec, 'NEGA', krr, ptstep, ppabst, ptht, prt, pths, prs, psvs, prhodj )
!
!
!-------------------------------------------------------------------------------
!
!*       4.     CLOUD ELECTRICITY
!               -----------------
!
!++cb++ 01/06/23
!IF (HELEC == 'ELE4') &
IF (HELEC(1:3) == 'ELE') THEN 
!--cb--
!
!*       4.1    Ion source from drift motion and cosmic rays
!
  CALL ION_SOURCE_ELEC (KTCOUNT, KRR, HLBCX, HLBCY,          &
                        PRHODREF, PRHODJ, PRT,               &
                        PSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND), &
                        PSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND), &
                        XEFIELDU, XEFIELDV, XEFIELDW         )
!
!*       4.2    Compute the coefficient that modifies the efficiency of IAGGS
!
  ALLOCATE(ZLATHAM_IAGGS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
  IF (LIAGGS_LATHAM) THEN
    ZLATHAM_IAGGS(:,:,:) = 1.0 + 0.4E-10 * MIN( 2.25E10,                &
                           XEFIELDU(:,:,:)**2+XEFIELDV(:,:,:)**2+XEFIELDW(:,:,:)**2 )
  ELSE
    ZLATHAM_IAGGS(:,:,:) = 1.0
  END IF
ELSE
  ALLOCATE(ZLATHAM_IAGGS(0,0,0))
END IF
!
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HCLOUD )
  CASE ('REVE')
!
!*       4.     REVERSIBLE MICROPHYSICAL SCHEME
!               -------------------------------
!
    CALL FAST_TERMS ( KRR, KMI, HRAD, HTURBDIM,                                &
                      HSCONV, HMF_CLOUD, OSUBG_COND, PTSTEP,                   &
                      PRHODJ, PSIGS, PPABST,                                   &
                      PCF_MF,PRC_MF,                                           &
                      PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),                    &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR                    )
!
  CASE ('KESS')
!
!*       5.     KESSLER MICROPHYSICAL SCHEME
!               ----------------------------
!
!
!*       5.1    Compute the explicit microphysical sources
!
    CALL SLOW_TERMS ( KSPLITR, PTSTEP, KMI, HSUBG_AUCV,                       &
                      PZZ, PRHODJ, PRHODREF, PCLDFR,                          &
                      PTHT, PRT(:,:,:,1), PRT(:,:,:,2), PRT(:,:,:,3), PPABST, &
                      PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),         &
                      PINPRR, PINPRR3D, PEVAP3D                         )
!
!*       5.2    Perform the saturation adjustment
!
    CALL FAST_TERMS ( KRR, KMI, HRAD, HTURBDIM,                                &
                      HSCONV, HMF_CLOUD, OSUBG_COND, PTSTEP,                   &
                      PRHODJ, PSIGS, PPABST,                                   &
                      PCF_MF,PRC_MF,                                           &
                      PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),                    &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2), PRRS=PRS(:,:,:,3), &
                      PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR                    )
!
!
  CASE ('C2R2','KHKO')
!
!*       7.     2-MOMENT WARM MICROPHYSICAL SCHEME C2R2 or KHKO
!               ---------------------------------------
!
!
!*       7.1    Compute the explicit microphysical sources
!
!
    CALL RAIN_C2R2_KHKO ( HCLOUD, OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP, KMI,     &
                     TPFILE, PZZ, PRHODJ, PRHODREF, PEXNREF,                      &
                     PPABST, PTHT, PRT(:,:,:,1), PRT(:,:,:,2),  PRT(:,:,:,3),     &
                     PTHM, PRCM, PPABSTT,                                          &
                     PW_ACT,PDTHRAD,PTHS, PRS(:,:,:,1),PRS(:,:,:,2),PRS(:,:,:,3), &
                     PSVT(:,:,:,NSV_C2R2BEG),   PSVT(:,:,:,NSV_C2R2BEG+1),        &
                     PSVT(:,:,:,NSV_C2R2BEG+2), PSVS(:,:,:,NSV_C2R2BEG),          &
                     PSVS(:,:,:,NSV_C2R2BEG+1), PSVS(:,:,:,NSV_C2R2BEG+2),        &
                     PINPRC, PINPRR, PINPRR3D, PEVAP3D ,                          &
                     PSVT(:,:,:,:), PSOLORG, PMI, HACTCCN,                        &
                     PINDEP, PSUPSAT, PNACT                                       )
!
!
!*       7.2    Perform the saturation adjustment
!
   IF (LSUPSAT) THEN
    CALL KHKO_NOTADJUST (KRR, KTCOUNT,TPFILE, HRAD,                              &
                         PTSTEP, PRHODJ, PPABSTT, PPABST, PRHODREF, PZZ,          &
                         PTHT,PRT(:,:,:,1),PRT(:,:,:,2),PRT(:,:,:,3),            &
                         PTHS,PRS(:,:,:,1),PRS(:,:,:,2),PRS(:,:,:,3),            &
                         PSVS(:,:,:,NSV_C2R2BEG+1), PSVS(:,:,:,NSV_C2R2BEG),     &
                         PSVS(:,:,:,NSV_C2R2BEG+3), PCLDFR, PSRCS, PNPRO, PSSPRO )
!
   ELSE
    CALL C2R2_ADJUST ( KRR,TPFILE, HRAD,                                &
                       HTURBDIM, OSUBG_COND, PTSTEP,                    &
                       PRHODJ, PSIGS, PPABST,                           &
                       PTHS=PTHS, PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2), &
                       PCNUCS=PSVS(:,:,:,NSV_C2R2BEG),                  &
                       PCCS=PSVS(:,:,:,NSV_C2R2BEG+1),                  &
                       PSRCS=PSRCS, PCLDFR=PCLDFR, PRRS=PRS(:,:,:,3)    )
!
   END IF
!
  CASE ('ICE3')
!
!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!$acc kernels
!$mnh_expand_array(JII=1:JIU,JJI=1:JJU,JKI=1:JKU)
    ZEXN(:,:,:)= (PPABST(:,:,:)/CST%XP00)**(CST%XRD/CST%XCPD)
!$mnh_end_expand_array(JII=1:JIU,JJI=1:JJU,JKI=1:JKU)
    IF (HELEC == 'ELE4') THEN
      ALLOCATE( ZCND    (SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)) )
      ALLOCATE( ZDEP    (SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)) )
      ALLOCATE( ZRCS_BEF(SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)) )
      ALLOCATE( ZRIS_BEF(SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)) )
    END IF
!
!*       9.1    Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZDZZ(:,:,1) = ZDZZ(:,:,IKB)
    ZDZZ(:,:,IKE+1) = ZDZZ(:,:,IKE)
!$acc end kernels
#ifndef MNH_OPENACC
    ZZZ(:,:,:) = MZF( PZZ(:,:,:) )
#else
    CALL MZF_DEVICE( PZZ, ZZZ )
#endif
    IF(LRED .AND. LADJ_BEFORE) THEN
      IF (HELEC == 'ELE4') THEN
        ! save the cloud droplets and ice crystals m.r. source before adjustement
        ZRCS_BEF(:,:,:) = PRS(:,:,:,2)
        ZRIS_BEF(:,:,:) = PRS(:,:,:,4)
      END IF
      !
      ! Performe the saturation ajdustment
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,           &
                      PARAM_ICEN, TBUCONF, KRR,                                &
                      'ADJU',                                                  &
                      PTSTEP, ZSIGQSAT2D,                                      &
                      PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV,PMFCONV, PPABST, ZZZ,  &
                      ZEXN, PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,          &
                      ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                            &
                      PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,        &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTH=PTHS*PTSTEP, PTHS=PTHS,                              &
                      OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR,    &
                      PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                      PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                      PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                      PRG=PRS(:,:,:,6)*PTSTEP,                                 &
                      TBUDGETS=TBUDGETS_PTR,KBUDGETS=SIZE(TBUDGETS_PTR), &
                      PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                    &
                      PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,                    &
                      PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF,       &
                      PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF)
      !
      IF (HELEC == 'ELE4') THEN
        ! Compute the condensation and sublimation rates
        ZCND(:,:,:) = PRS(:,:,:,2) - ZRCS_BEF(:,:,:)
        ZDEP(:,:,:) = PRS(:,:,:,4) - ZRIS_BEF(:,:,:)
        !
        ! Compute the charge exchanged during evaporation of cloud droplets (negative ZCND) and
        !                              during sublimation of ice crystals (negative ZDEP)
        CALL ELEC_ADJUST (KRR, PRHODJ, HCLOUD, 'ADJU',                                   &
                          PRC=ZRCS_BEF(:,:,:)*PTSTEP, PRI=ZRIS_BEF(:,:,:)*PTSTEP,        &
                          PQC=PSVS(:,:,:,NSV_ELECBEG+1)*PTSTEP,                          &
                          PQI=PSVS(:,:,:,NSV_ELECBEG+3)*PTSTEP,                          &
                          PQCS=PSVS(:,:,:,NSV_ELECBEG+1), PQIS=PSVS(:,:,:,NSV_ELECBEG+3),&
                          PQPIS=PSVS(:,:,:,NSV_ELECBEG), PQNIS=PSVS(:,:,:,NSV_ELECEND),  &
                          PCND=ZCND, PDEP=ZDEP)
      END IF
    ENDIF
    IF (LRED) THEN
      IF (HELEC == 'ELE4') THEN
        ! to match with PHYEX, electric charge variables are no more optional, but their size
        ! depends on the activation (or not) of the electrification scheme
        GELEC = .TRUE.
        ALLOCATE(ZQPIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQNIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQCT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQRT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQST(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQGT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQPIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQNIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQCS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQRS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQSS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQGS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ZQPIT(:,:,:) = PSVT(:,:,:,NSV_ELECBEG)
        ZQCT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+1)
        ZQRT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+2)
        ZQIT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+3)
        ZQST(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+4)
        ZQGT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+5)
        ZQNIT(:,:,:) = PSVT(:,:,:,NSV_ELECEND)
        ZQPIS(:,:,:) = PSVS(:,:,:,NSV_ELECBEG)
        ZQCS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+1)
        ZQRS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+2)
        ZQIS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+3)
        ZQSS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+4)
        ZQGS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+5)
        ZQNIS(:,:,:) = PSVS(:,:,:,NSV_ELECEND)
        IF (ELEC_DESCR%LSEDIM_BEARD) THEN
          ALLOCATE(ZEFIELDW(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
          ZEFIELDW(:,:,:) = XEFIELDW(:,:,:)
        ELSE
          ALLOCATE(ZEFIELDW(0,0,0))
        END IF
      ELSE
        GELEC = .FALSE.
        ALLOCATE(ZQPIT(0,0,0))
        ALLOCATE(ZQNIT(0,0,0))
        ALLOCATE(ZQCT(0,0,0))
        ALLOCATE(ZQRT(0,0,0))
        ALLOCATE(ZQIT(0,0,0))
        ALLOCATE(ZQST(0,0,0))
        ALLOCATE(ZQGT(0,0,0))
        ALLOCATE(ZQPIS(0,0,0))
        ALLOCATE(ZQNIS(0,0,0))
        ALLOCATE(ZQCS(0,0,0))
        ALLOCATE(ZQRS(0,0,0))
        ALLOCATE(ZQIS(0,0,0))
        ALLOCATE(ZQSS(0,0,0))
        ALLOCATE(ZQGS(0,0,0))
        ALLOCATE(ZEFIELDW(0,0,0))
      END IF
      ALLOCATE(ZQHT(0,0,0))
      ALLOCATE(ZQHS(0,0,0))
      !
      CALL RAIN_ICE (YLDIMPHYEX,CST, PARAM_ICEN, RAIN_ICE_PARAMN, RAIN_ICE_DESCRN, &
                    ELEC_PARAM, ELEC_DESCR, TBUCONF, GELEC, ELEC_DESCR%LSEDIM_BEARD,       &
                    XTHVREFZ(IKB),                                              &
                    PTSTEP, KRR, ZEXN,                                          &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,      &
                    ZDUM, ZDUM, ZDUM, ZDUM, &
                    PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                     &
                    PTHT, PRT, PTHS, PRS, &
                    PINPRC,PINPRR, PEVAP3D,                                     &
                    PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,                     &
                    TBUDGETS_PTR, SIZE(TBUDGETS_PTR), &
                    ZQPIT, ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQNIT,                 &
                    ZQPIS, ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQNIS,                 &
                    ZEFIELDW, ZLATHAM_IAGGS,                                    &
                    PSEA,PTOWN, PFPR=ZFPR                                       )
      !
      IF (HELEC == 'ELE4') THEN
        PSVT(:,:,:,NSV_ELECBEG)   = ZQPIT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+1) = ZQCT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+2) = ZQRT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+3) = ZQIT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+4) = ZQST(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+5) = ZQGT(:,:,:)
        PSVT(:,:,:,NSV_ELECEND)   = ZQNIT(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG)   = ZQPIS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+1) = ZQCS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+2) = ZQRS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+3) = ZQIS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+4) = ZQSS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+5) = ZQGS(:,:,:)
        PSVS(:,:,:,NSV_ELECEND)   = ZQNIS(:,:,:)
      END IF
      DEALLOCATE(ZQPIT)
      DEALLOCATE(ZQNIT)
      DEALLOCATE(ZQCT)
      DEALLOCATE(ZQRT)
      DEALLOCATE(ZQIT)
      DEALLOCATE(ZQST)
      DEALLOCATE(ZQGT)
      DEALLOCATE(ZQHT)
      DEALLOCATE(ZQPIS)
      DEALLOCATE(ZQNIS)
      DEALLOCATE(ZQCS)
      DEALLOCATE(ZQRS)
      DEALLOCATE(ZQIS)
      DEALLOCATE(ZQSS)
      DEALLOCATE(ZQGS)
      DEALLOCATE(ZQHS)
      DEALLOCATE(ZEFIELDW)
      !
    ELSE 
      IF (HELEC == 'ELE3') THEN
        ! --> old version of the electrification scheme
        ! Should be removed in a future version of MNH once the new electrification scheme is fully validated
        ! Compute the explicit microphysical sources and the explicit charging rates
        CALL RAIN_ICE_ELEC (OSEDIC, HSUBG_AUCV, OWARM,                            &
                            KSPLITR, PTSTEP, KMI, KRR,                            &
                            PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                            PTHT, PRT(:,:,:,1), PRT(:,:,:,2), PRT(:,:,:,3),       &
                            PRT(:,:,:,4), PRT(:,:,:,5), PRT(:,:,:,6),             &
                            PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                            PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                            PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                            PINPRS, PINPRG, PSIGS,                                &
                            PSVT(:,:,:,NSV_ELECBEG),   PSVT(:,:,:,NSV_ELECBEG+1), &
                            PSVT(:,:,:,NSV_ELECBEG+2), PSVT(:,:,:,NSV_ELECBEG+3), &
                            PSVT(:,:,:,NSV_ELECBEG+4), PSVT(:,:,:,NSV_ELECBEG+5), &
                            PSVT(:,:,:,NSV_ELECEND),                              &
                            PSVS(:,:,:,NSV_ELECBEG),   PSVS(:,:,:,NSV_ELECBEG+1), &
                            PSVS(:,:,:,NSV_ELECBEG+2), PSVS(:,:,:,NSV_ELECBEG+3), &
                            PSVS(:,:,:,NSV_ELECBEG+4), PSVS(:,:,:,NSV_ELECBEG+5), &
                            PSVS(:,:,:,NSV_ELECEND),                              &
                            PSEA, PTOWN                                           )
      ELSE
        CALL RAIN_ICE_OLD (YLDIMPHYEX, OSEDIC, CSEDIM, HSUBG_AUCV, OWARM, 1, IKU, 1,&
                      KSPLITR, PTSTEP, KRR,                                         &
                      ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,        &
                      PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                             &
                      PRT(:,:,:,3), PRT(:,:,:,4),                                   &
                      PRT(:,:,:,5), PRT(:,:,:,6),                                   &
                      PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),               &
                      PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),                     &
                      PINPRC,PINPRR, PINPRR3D, PEVAP3D,                             &
                      PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,                        &
                      PSEA, PTOWN, PFPR=ZFPR                                        )
      END IF
    END IF

!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
      IF (HELEC == 'ELE4') THEN
        ! save the cloud droplets and ice crystals m.r. source before adjustement
        ZRCS_BEF(:,:,:) = PRS(:,:,:,2)
        ZRIS_BEF(:,:,:) = PRS(:,:,:,4)
      END IF
      !
      ! Perform the saturation ajdustment
      IF (HELEC == 'ELE3') THEN
        ! --> old version of the electrification scheme
        CALL ICE_ADJUST_ELEC (KRR, KMI, HRAD, HTURBDIM,                               &
                              HSCONV, HMF_CLOUD,                                      &
                              OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,                   &
                              PRHODJ, PEXNREF, PSIGS, PPABST, ZZZ,                    &
                              PMFCONV, PCF_MF, PRC_MF, PRI_MF,                        &
                              PRT(:,:,:,1), PRT(:,:,:,2), PRS(:,:,:,1), PRS(:,:,:,2), &
                              PTHS, PSRCS, PCLDFR,                                    &
                              PRT(:,:,:,3), PRS(:,:,:,3), PRT(:,:,:,4), PRS(:,:,:,4), &
                              PRT(:,:,:,5), PRS(:,:,:,5), PRT(:,:,:,6), PRS(:,:,:,6), &
                              PSVT(:,:,:,NSV_ELECBEG),   PSVS(:,:,:,NSV_ELECBEG),     &
                              PSVT(:,:,:,NSV_ELECBEG+1), PSVS(:,:,:,NSV_ELECBEG+1),   &
                              PSVT(:,:,:,NSV_ELECBEG+2), PSVS(:,:,:,NSV_ELECBEG+2),   &
                              PSVT(:,:,:,NSV_ELECBEG+3), PSVS(:,:,:,NSV_ELECBEG+3),   &
                              PSVT(:,:,:,NSV_ELECBEG+4), PSVS(:,:,:,NSV_ELECBEG+4),   &
                              PSVT(:,:,:,NSV_ELECBEG+5), PSVS(:,:,:,NSV_ELECBEG+5),   &
                              PSVT(:,:,:,NSV_ELECEND),   PSVS(:,:,:,NSV_ELECEND)      )
      ELSE
        CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,                   & 
                         PARAM_ICEN, TBUCONF, KRR, 'DEPI',                               &
                         PTSTEP, ZSIGQSAT2D,                                             &
                         PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,PPABST, ZZZ, &
                         ZEXN, PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,                 &
                         ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                                   &
                         PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,               &
                         PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                           &
                         PTH=PTHS*PTSTEP, PTHS=PTHS,                                     &
                         OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR,     &
                         PRR=PRS(:,:,:,3)*PTSTEP,                                        &
                         PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),                     &
                         PRS=PRS(:,:,:,5)*PTSTEP,                                        &
                         PRG=PRS(:,:,:,6)*PTSTEP,                                        &
                         TBUDGETS=TBUDGETS_PTR, KBUDGETS=SIZE(TBUDGETS_PTR), &
                         PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                           &
                         PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,                           &
                         PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF,               &
                         PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF                )

        !
        IF (HELEC == 'ELE4') THEN
          ! Compute the condensation and sublimation rates
          ZCND(:,:,:) = PRS(:,:,:,2) - ZRCS_BEF(:,:,:)
          ZDEP(:,:,:) = PRS(:,:,:,4) - ZRIS_BEF(:,:,:)
          !
          ! Compute the charge exchanged during evaporation of cloud droplets (negative ZCND) and
          !                              during sublimation of ice crystals (negative ZDEP)
          CALL ELEC_ADJUST (KRR, PRHODJ, HCLOUD, 'DEPI',                                   &
                            PRC=ZRCS_BEF(:,:,:)*PTSTEP, PRI=ZRIS_BEF(:,:,:)*PTSTEP,        &
                            PQC=PSVS(:,:,:,NSV_ELECBEG+1)*PTSTEP,                          &
                            PQI=PSVS(:,:,:,NSV_ELECBEG+3)*PTSTEP,                          &
                            PQCS=PSVS(:,:,:,NSV_ELECBEG+1), PQIS=PSVS(:,:,:,NSV_ELECBEG+3),&
                            PQPIS=PSVS(:,:,:,NSV_ELECBEG), PQNIS=PSVS(:,:,:,NSV_ELECEND),  &
                            PCND=ZCND, PDEP=ZDEP)
        END IF
      END IF
    END IF
!
  CASE ('ICE4')
!
!*       10.    MIXED-PHASE MICROPHYSICAL SCHEME (WITH 4 ICE SPECIES)
!               -----------------------------------------------------
!
!$acc kernels
 !$mnh_expand_array(JII=1:JIU,JJI=1:JJU,JKI=1:JKU)
    ZEXN(:,:,:)= (PPABST(:,:,:)/CST%XP00)**(CST%XRD/CST%XCPD)
 !$mnh_end_expand_array(JII=1:JIU,JJI=1:JJU,JKI=1:JKU)
!
!*       10.1   Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
!$acc end kernels
#ifndef MNH_OPENACC
    ZZZ(:,:,:) = MZF( PZZ(:,:,:) )
#else
    CALL MZF_DEVICE( PZZ, ZZZ )
#endif
    IF(LRED .AND. LADJ_BEFORE) THEN
      IF (HELEC == 'ELE4') THEN
        ! save the cloud droplets and ice crystals m.r. source before adjustement
        ZRCS_BEF(:,:,:) = PRS(:,:,:,2)
        ZRIS_BEF(:,:,:) = PRS(:,:,:,4)
      END IF
      !
      ! Perform the saturation ajdustment
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,           &
                       PARAM_ICEN, TBUCONF, KRR,                               &
                       'ADJU',                                                 &
                       PTSTEP, ZSIGQSAT2D,                                     &
                       PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV,PMFCONV, PPABST, ZZZ, &
                       ZEXN, PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,         &
                       ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                           &
                       PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,       &
                       PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                   &
                       PTH=PTHS*PTSTEP, PTHS=PTHS,                             &
                       OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR, &
                       PRR=PRS(:,:,:,3)*PTSTEP,                                &
                       PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),             &
                       PRS=PRS(:,:,:,5)*PTSTEP,                                &
                       PRG=PRS(:,:,:,6)*PTSTEP,                                &
                       TBUDGETS=TBUDGETS_PTR, KBUDGETS=SIZE(TBUDGETS_PTR), &
                       PRH=PRS(:,:,:,7)*PTSTEP,                                &
                       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                   &
                       PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,                   &
                       PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF,       &
                       PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF        )
      !
      IF (HELEC == 'ELE4') THEN
        ! Compute the condensation and sublimation rates
        ZCND(:,:,:) = PRS(:,:,:,2) - ZRCS_BEF(:,:,:)
        ZDEP(:,:,:) = PRS(:,:,:,4) - ZRIS_BEF(:,:,:)
        !
        ! Compute the charge exchanged during evaporation of cloud droplets (negative ZCND) and
        !                              during sublimation of ice crystals (negative ZDEP)
        CALL ELEC_ADJUST (KRR, PRHODJ, HCLOUD, 'ADJU',                                   &
                          PRC=ZRCS_BEF(:,:,:)*PTSTEP, PRI=ZRIS_BEF(:,:,:)*PTSTEP,        &
                          PQC=PSVS(:,:,:,NSV_ELECBEG+1)*PTSTEP,                          &
                          PQI=PSVS(:,:,:,NSV_ELECBEG+3)*PTSTEP,                          &
                          PQCS=PSVS(:,:,:,NSV_ELECBEG+1), PQIS=PSVS(:,:,:,NSV_ELECBEG+3),&
                          PQPIS=PSVS(:,:,:,NSV_ELECBEG), PQNIS=PSVS(:,:,:,NSV_ELECEND),  &
                          PCND=ZCND, PDEP=ZDEP                                           )
      END IF
    ENDIF
    IF  (LRED) THEN
      IF (HELEC == 'ELE4') THEN
        ! to match with PHYEX, electric charge variables are no more optional, but their size
        ! depends on the activation (or not) of the electrification scheme
        GELEC = .TRUE.
        ALLOCATE(ZQPIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQNIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQCT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQRT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQIT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQST(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQGT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQHT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQPIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQNIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQCS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQRS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQIS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQSS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQGS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ALLOCATE(ZQHS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
        ZQPIT(:,:,:) = PSVT(:,:,:,NSV_ELECBEG)
        ZQCT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+1)
        ZQRT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+2)
        ZQIT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+3)
        ZQST(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+4)
        ZQGT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+5)
        ZQHT(:,:,:)  = PSVT(:,:,:,NSV_ELECBEG+6)
        ZQNIT(:,:,:) = PSVT(:,:,:,NSV_ELECEND)
        ZQPIS(:,:,:) = PSVS(:,:,:,NSV_ELECBEG)
        ZQCS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+1)
        ZQRS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+2)
        ZQIS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+3)
        ZQSS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+4)
        ZQGS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+5)
        ZQHS(:,:,:)  = PSVS(:,:,:,NSV_ELECBEG+6)
        ZQNIS(:,:,:) = PSVS(:,:,:,NSV_ELECEND)
        IF (ELEC_DESCR%LSEDIM_BEARD) THEN
          ALLOCATE(ZEFIELDW(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
          ZEFIELDW(:,:,:) = XEFIELDW(:,:,:)
        ELSE
          ALLOCATE(ZEFIELDW(0,0,0))
        END IF
      ELSE
        GELEC = .FALSE.
        ALLOCATE(ZQPIT(0,0,0))
        ALLOCATE(ZQNIT(0,0,0))
        ALLOCATE(ZQCT(0,0,0))
        ALLOCATE(ZQRT(0,0,0))
        ALLOCATE(ZQIT(0,0,0))
        ALLOCATE(ZQST(0,0,0))
        ALLOCATE(ZQGT(0,0,0))
        ALLOCATE(ZQHT(0,0,0))
        ALLOCATE(ZQPIS(0,0,0))
        ALLOCATE(ZQNIS(0,0,0))
        ALLOCATE(ZQCS(0,0,0))
        ALLOCATE(ZQRS(0,0,0))
        ALLOCATE(ZQIS(0,0,0))
        ALLOCATE(ZQSS(0,0,0))
        ALLOCATE(ZQGS(0,0,0))
        ALLOCATE(ZQHS(0,0,0))
        ALLOCATE(ZEFIELDW(0,0,0))
      END IF
      !
      CALL RAIN_ICE (YLDIMPHYEX,CST, PARAM_ICEN, RAIN_ICE_PARAMN, RAIN_ICE_DESCRN, &
                     ELEC_PARAM, ELEC_DESCR, TBUCONF, GELEC, ELEC_DESCR%LSEDIM_BEARD,        &
                     XTHVREFZ(IKB),                                              &
                     PTSTEP, KRR, ZEXN,                                          &
                     ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,      &
                     ZDUM, ZDUM, ZDUM, ZDUM, &
                     PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                     &
                     PTHT, PRT, PTHS, PRS, &
                     PINPRC, PINPRR, PEVAP3D,                                    &
                     PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,                     &
                     TBUDGETS_PTR, SIZE(TBUDGETS_PTR), &
                     ZQPIT, ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQNIT,                 &
                     ZQPIS, ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQNIS,                 &
                     ZEFIELDW, ZLATHAM_IAGGS,                                    &
                     PSEA, PTOWN,                                                &
                     PINPRH, PFPR=ZFPR, &
                     PQHT=ZQHT, PQHS=ZQHS                                        )
      !
      IF (HELEC == 'ELE4') THEN
        PSVT(:,:,:,NSV_ELECBEG)   = ZQPIT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+1) = ZQCT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+2) = ZQRT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+3) = ZQIT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+4) = ZQST(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+5) = ZQGT(:,:,:)
        PSVT(:,:,:,NSV_ELECBEG+6) = ZQHT(:,:,:)
        PSVT(:,:,:,NSV_ELECEND)   = ZQNIT(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG)   = ZQPIS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+1) = ZQCS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+2) = ZQRS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+3) = ZQIS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+4) = ZQSS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+5) = ZQGS(:,:,:)
        PSVS(:,:,:,NSV_ELECBEG+6) = ZQHS(:,:,:)
        PSVS(:,:,:,NSV_ELECEND)   = ZQNIS(:,:,:)
      END IF
      DEALLOCATE(ZQPIT)
      DEALLOCATE(ZQNIT)
      DEALLOCATE(ZQCT)
      DEALLOCATE(ZQRT)
      DEALLOCATE(ZQIT)
      DEALLOCATE(ZQST)
      DEALLOCATE(ZQGT)
      DEALLOCATE(ZQHT)
      DEALLOCATE(ZQPIS)
      DEALLOCATE(ZQNIS)
      DEALLOCATE(ZQCS)
      DEALLOCATE(ZQRS)
      DEALLOCATE(ZQIS)
      DEALLOCATE(ZQSS)
      DEALLOCATE(ZQGS)
      DEALLOCATE(ZQHS)
      DEALLOCATE(ZEFIELDW)
      !
    ELSE
      CALL RAIN_ICE_OLD (YLDIMPHYEX, OSEDIC, CSEDIM, HSUBG_AUCV, OWARM, 1, IKU, 1,    &
                    KSPLITR, PTSTEP, KRR,                                 &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                    PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,                &
                    PSEA, PTOWN,                                          &
                    PRT(:,:,:,7), PRS(:,:,:,7), PINPRH, PFPR=ZFPR         )
    END IF
!
!
!*       10.2   Perform the saturation adjustment over cloud ice and cloud water
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
      IF (HELEC == 'ELE4') THEN
        ! save the cloud droplets and ice crystals m.r. source before adjustement
        ZRCS_BEF(:,:,:) = PRS(:,:,:,2)
        ZRIS_BEF(:,:,:) = PRS(:,:,:,4)
      END IF
      !
      ! Perform the saturation ajdustment
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,                  & 
                      PARAM_ICEN, TBUCONF, KRR, 'DEPI',                               &
                      PTSTEP, ZSIGQSAT2D,                                             &
                      PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,PPABST, ZZZ, &
                      ZEXN, PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,                 &
                      ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                                   &
                      PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,               &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                           &
                      PTH=PTHS*PTSTEP, PTHS=PTHS,                                     &
                      OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR,     &
                      PRR=PRS(:,:,:,3)*PTSTEP,                                        &
                      PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),                     &
                      PRS=PRS(:,:,:,5)*PTSTEP,                                        &
                      PRG=PRS(:,:,:,6)*PTSTEP,                                        &
                      TBUDGETS=TBUDGETS_PTR, KBUDGETS=SIZE(TBUDGETS_PTR), &
                      PRH=PRS(:,:,:,7)*PTSTEP,                                        &
                      PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                           &
                      PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,                           &
                      PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF,               &
                      PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF                )

      !
      IF (HELEC == 'ELE4') THEN
        ! Compute the condensation and sublimation rates
        ZCND(:,:,:) = PRS(:,:,:,2) - ZRCS_BEF(:,:,:)
        ZDEP(:,:,:) = PRS(:,:,:,4) - ZRIS_BEF(:,:,:)
        !
        ! Compute the charge exchanged during evaporation of cloud droplets (negative ZCND) and
        !                              during sublimation of ice crystals (negative ZDEP)
        CALL ELEC_ADJUST (KRR, PRHODJ, HCLOUD, 'ADJU',                                   &
                          PRC=ZRCS_BEF(:,:,:)*PTSTEP, PRI=ZRIS_BEF(:,:,:)*PTSTEP,        &
                          PQC=PSVS(:,:,:,NSV_ELECBEG+1)*PTSTEP,                          &
                          PQI=PSVS(:,:,:,NSV_ELECBEG+3)*PTSTEP,                          &
                          PQCS=PSVS(:,:,:,NSV_ELECBEG+1), PQIS=PSVS(:,:,:,NSV_ELECBEG+3),&
                          PQPIS=PSVS(:,:,:,NSV_ELECBEG), PQNIS=PSVS(:,:,:,NSV_ELECEND),  &
                          PCND=ZCND, PDEP=ZDEP                                           )
      END IF
    END IF

!           
!
!*       12.    2-MOMENT MIXED-PHASE MICROPHYSICAL SCHEME LIMA
!               --------------------------------------------------------------
!
!
!*       12.1   Compute the explicit microphysical sources
!
  CASE ('LIMA')
    !
    IF (HELEC == 'ELE4') THEN
      GELEC = .TRUE.
    ELSE
      GELEC = .FALSE.
    END IF
    !
    LLDTHRAD = SIZE(PDTHRAD) /= 0
!$acc kernels
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
!$acc end kernels
#ifndef MNH_OPENACC
    ZZZ(:,:,:) = MZF( PZZ(:,:,:) )
#else
    CALL MZF_DEVICE( PZZ, ZZZ )
#endif
    IF (LPTSPLIT) THEN 
      IF (GELEC) THEN
         CALL LIMA (PARAM_LIMA, PARAM_LIMA_WARM, PARAM_LIMA_COLD, PARAM_LIMA_MIXED,&
                   TNSV, YLDIMPHYEX,CST, NEBN, RAIN_ICE_DESCRN, RAIN_ICE_PARAMN,       &
                   ELEC_DESCR, ELEC_PARAM,                                 &
                   TBUCONF, TBUDGETS_PTR, HACTCCN, SIZE(TBUDGETS_PTR), KRR, &
                   PTSTEP, GELEC,                                          &
                   PRHODREF, PEXNREF, ZDZZ, XTHVREFZ(IKB),                 &
                   PRHODJ, PPABST,                                         &
                   NCARB, NSOA, NSP, LDUST, LSALT, LORILAM,                &
                   LLDTHRAD, PDTHRAD, PTHT, PRT,                           &
                   PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), PCIT, PW_ACT,    &
                   PSVT, PSOLORG, PMI,                                     &
                   PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),       &
                   PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, PINPRH, &
                   PEVAP3D, PCLDFR, PICEFR, PRAINFR, ZFPR,                 &
                   PHLC_HCF, PHLC_HRC,                                     &
                   PHLI_HCF, PHLI_HRI,                                     &
                   ZLATHAM_IAGGS, XEFIELDW,                                &
                   PSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),                    &
                   PSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND)                     )
      ELSE
        CALL LIMA (PARAM_LIMA, PARAM_LIMA_WARM, PARAM_LIMA_COLD, PARAM_LIMA_MIXED,&
                   TNSV, YLDIMPHYEX,CST, NEBN, RAIN_ICE_DESCRN, RAIN_ICE_PARAMN,       &
                   ELEC_DESCR, ELEC_PARAM,                                 &
                   TBUCONF, TBUDGETS_PTR, HACTCCN, SIZE(TBUDGETS_PTR), KRR, &
                   PTSTEP, GELEC,                                          &
                   PRHODREF, PEXNREF, ZDZZ, XTHVREFZ(IKB),                 &
                   PRHODJ, PPABST,                                         &
                   NCARB, NSOA, NSP, LDUST, LSALT, LORILAM,                &
                   LLDTHRAD, PDTHRAD, PTHT, PRT,                           &
                   PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), PCIT, PW_ACT,    &
                   PSVT, PSOLORG, PMI,                                     &
                   PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),       &
                   PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, PINPRH, &
                   PEVAP3D, PCLDFR, PICEFR, PRAINFR, ZFPR,                 &
                   PHLC_HCF, PHLC_HRC,                                     &
                   PHLI_HCF, PHLI_HRI,                                     &
                   ZLATHAM_IAGGS                                           )
      END IF
    ELSE
      IF (OWARM) CALL LIMA_WARM(OACTIT, HACTCCN, OSEDC, ORAIN, KSPLITR, PTSTEP,   &
                                KMI, TPFILE, KRR, PZZ, PRHODJ,                    &
                                PRHODREF, PEXNREF, PW_ACT, PPABST,                &
                                PDTHRAD,                                          &
                                PTHT, PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                PSVT, PSOLORG, PMI,                               &
                                PINPRC, PINPRR, PINDEP, PINPRR3D, PEVAP3D         )
!
      IF (NMOM_I.GE.1) CALL LIMA_COLD(CST, OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,         &
                                      KRR, PZZ, PRHODJ,                                 &
                                      PRHODREF, PEXNREF, PPABST, PW_ACT,                &
                                      PTHT, PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                      PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                      PINPRS, PINPRG, PINPRH                            )
!
      IF (OWARM .AND. NMOM_I.GE.1) CALL LIMA_MIXED(OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,              &
                                                   KRR, PZZ, PRHODJ,                                 &
                                                   PRHODREF, PEXNREF, PPABST, PW_ACT,                &
                                                   PTHT, PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                                   PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END)  )
    ENDIF
!
!*       12.2   Perform the saturation adjustment
!
    IF (HELEC == 'ELE4') THEN
      ! save the cloud droplets and ice crystals m.r. source before adjustement
      ZRCS_BEF(:,:,:) = PRS(:,:,:,2)
      ZRIS_BEF(:,:,:) = PRS(:,:,:,4)
    END IF
    !
    IF (LSPRO) THEN
      CALL LIMA_NOTADJUST (KMI, TPFILE, HRAD,                                       &
                           PTSTEP, PRHODJ, PPABSTT, PPABST, PRHODREF, PEXNREF, PZZ, &
                           PTHT,PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),         &
                           PTHS,PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),         &
                           PCLDFR, PICEFR, PRAINFR, PSRCS                           )
   ELSE IF (LPTSPLIT) THEN
      ! currently using ZSIGQSAT2D as a dummy argument for PICE_CLD_WGT, only used in condensation with OCND2
       CALL LIMA_ADJUST_SPLIT(PARAM_LIMA, PARAM_LIMA_WARM, &
                             TNSV, YLDIMPHYEX, CST, NEBN, TURBN, TBUCONF, TBUDGETS_PTR, SIZE(TBUDGETS_PTR), &
                             KRR, CCONDENS, CLAMBDA3,                                        &
                             NCARB, NSOA, NSP, LDUST, LSALT, LORILAM,                        &
                             OSUBG_COND, OSIGMAS, PTSTEP, ZSIGQSAT2D,                        &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS,                               &
                             SIZE(PMFCONV)/=0, PMFCONV, PPABST, ZZZ,                         &
                             LLDTHRAD, PDTHRAD, PW_ACT,                                      &
                             PRT, PRS, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),                &
                             PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),                          &
                             HACTCCN, PSVT, PSOLORG, PMI,                                    &
                             PTHS, SIZE(PSRCS, 3)/=0, PSRCS, PCLDFR, PICEFR,                 &
                             PRC_MF, PRI_MF, PCF_MF,                                         &
                             ZSIGQSAT2D, PWEIGHT_MF_CLOUD,                        &
                             PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                         &
                             PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF              )
    ELSE
      CALL LIMA_ADJUST(KRR, KMI, TPFILE,                                &
                       OSUBG_COND, PTSTEP,                              &
                       PRHODREF, PRHODJ, PEXNREF, PPABST, PPABSTT,      &
                       PRT, PRS, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                       PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),           &
                       PTHS, PSRCS, PCLDFR, PICEFR, PRAINFR             )
    ENDIF
    !
    IF (HELEC == 'ELE4') THEN
      ! Compute the condensation and sublimation rates
      ZCND(:,:,:) = PRS(:,:,:,2) - ZRCS_BEF(:,:,:)
      ZDEP(:,:,:) = PRS(:,:,:,4) - ZRIS_BEF(:,:,:)
      ! Compute the charge exchanged during evaporation of cloud droplets (negative ZCND) and
      !                              during sublimation of ice crystals (negative ZDEP)
      CALL ELEC_ADJUST (KRR, PRHODJ, HCLOUD, 'CEDS',                                   &
                        PRC=ZRCS_BEF(:,:,:)*PTSTEP, PRI=ZRIS_BEF(:,:,:)*PTSTEP,        &
                        PQC=PSVS(:,:,:,NSV_ELECBEG+1)*PTSTEP,                          &
                        PQI=PSVS(:,:,:,NSV_ELECBEG+3)*PTSTEP,                          &
                        PQCS=PSVS(:,:,:,NSV_ELECBEG+1), PQIS=PSVS(:,:,:,NSV_ELECBEG+3),&
                        PQPIS=PSVS(:,:,:,NSV_ELECBEG), PQNIS=PSVS(:,:,:,NSV_ELECEND),  &
                        PCND=ZCND, PDEP=ZDEP                                           )
    END IF
!
END SELECT
!
IF (ALLOCATED(ZLATHAM_IAGGS)) DEALLOCATE(ZLATHAM_IAGGS)
!
IF(HCLOUD=='ICE3' .OR. HCLOUD=='ICE4' ) THEN
! TODO: code a generic routine to update vertical lower and upper levels to 0, a
! specific value or to IKB or IKE and apply it to every output prognostic variable of physics
!$acc kernels
  PCIT(:,:,1)     = 0.
  PCIT(:,:,IKE+1) = 0.

  !$mnh_expand_where(JII=1:JIU,JJI=1:JJU,JKI=1:JKU) 
  PINPRC3D(:,:,:)=ZFPR(:,:,:,2) / CST%XRHOLW
  PINPRR3D(:,:,:)=ZFPR(:,:,:,3) / CST%XRHOLW
  PINPRS3D(:,:,:)=ZFPR(:,:,:,5) / CST%XRHOLW
  PINPRG3D(:,:,:)=ZFPR(:,:,:,6) / CST%XRHOLW
  IF(KRR==7) THEN
    PINPRH3D(:,:,:)=ZFPR(:,:,:,7) / CST%XRHOLW
  END IF
  WHERE (PRT(:,:,:,2) > 1.E-04 )
    PSPEEDC(:,:,:)=ZFPR(:,:,:,2) / (PRT(:,:,:,2) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,3) > 1.E-04 )
    PSPEEDR(:,:,:)=ZFPR(:,:,:,3) / (PRT(:,:,:,3) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,5) > 1.E-04 )
    PSPEEDS(:,:,:)=ZFPR(:,:,:,5) / (PRT(:,:,:,5) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,6) > 1.E-04 )
    PSPEEDG(:,:,:)=ZFPR(:,:,:,6) / (PRT(:,:,:,6) * PRHODREF(:,:,:))
  ENDWHERE
  IF(KRR==7) THEN
    WHERE (PRT(:,:,:,7) > 1.E-04 )
      PSPEEDH(:,:,:)=ZFPR(:,:,:,7) / (PRT(:,:,:,7) * PRHODREF(:,:,:))
    ENDWHERE
  ENDIF
  !$mnh_end_expand_where(JII=1:JIU,JJI=1:JJU,JKI=1:JKU)
!$acc end kernels
ENDIF
!
! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
call Sources_neg_correct( hcloud, helec, 'NECON', krr, ptstep, ppabst, ptht, prt, pths, prs, psvs, prhodj )
!
!-------------------------------------------------------------------------------
!
!*      13.     SWITCH BACK TO THE PROGNOSTIC VARIABLES
!               ---------------------------------------
!
!$acc kernels
PTHS(:,:,:) = PTHS(:,:,:) * PRHODJ(:,:,:)
!
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) * PRHODJ(:,:,:)
END DO
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
  DO JSV = ISVBEG, ISVEND
    PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) * PRHODJ(:,:,:)
  ENDDO
ENDIF
!
IF (HELEC /= 'NONE') THEN
  DO JSV = NSV_ELECBEG, NSV_ELECEND
    PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) * PRHODJ(:,:,:)
  END DO
!
!++cb-- ce qui suit n'est plus present en version standard en 5-6 : pourquoi ?
! Note that the LiNOx Conc. (in mol/mol) is PSVS (:,::,NSV_LNOXBEG)
! but there is no need to *PRHODJ(:,:,:) as it is done implicitly
! during unit conversion in flash_geom.
!
  PSVS(:,:,:,NSV_ELECBEG) = MAX(0., PSVS(:,:,:,NSV_ELECBEG))
  PSVS(:,:,:,NSV_ELECEND) = MAX(0., PSVS(:,:,:,NSV_ELECEND))
END IF
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RESOLVED_CLOUD
