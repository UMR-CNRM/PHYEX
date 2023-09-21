!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_RESOLVED_CLOUD
!     ##########################
INTERFACE
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HACTCCN, HSCONV, HMF_CLOUD,                  &
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
                                  PCF_MF,PRC_MF, PRI_MF,                               &
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
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE   ! Output file
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
!     ##########################################################################
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HACTCCN, HSCONV, HMF_CLOUD,                  &
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
                                  PCF_MF,PRC_MF, PRI_MF,                               &
                                  PINPRC,PINPRC3D,PINPRR,PINPRR3D, PEVAP3D,            &
                                  PINPRS,PINPRS3D,PINPRG,PINPRG3D,PINPRH,PINPRH3D,     &
                                  PSOLORG,PMI,                                         &
                                  PSPEEDC, PSPEEDR, PSPEEDS, PSPEEDG, PSPEEDH,         &
                                  PINDEP, PSUPSAT,  PNACT, PNPRO,PSSPRO, PRAINFR,      &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,              &
                                  PSEA,PTOWN          )   
!     ##########################################################################
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_BUDGET,         ONLY: TBUDGETS, TBUCONF
USE MODD_CH_AEROSOL,     ONLY: LORILAM
USE MODD_DUST,           ONLY: LDUST
USE MODD_CST,            ONLY: CST
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_DUST ,          ONLY: LDUST
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_NEB_n,          ONLY: NEBN, CCONDENS, CLAMBDA3
USE MODD_NSV,            ONLY: NSV, NSV_C1R3END, NSV_C2R2BEG, NSV_C2R2END,                       &
                               NSV_LIMA_BEG, NSV_LIMA_END, NSV_LIMA_CCN_FREE, NSV_LIMA_IFN_FREE, &
                               NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_NR, NSV_AEREND,NSV_DSTEND,NSV_SLTEND
USE MODD_PARAM_C2R2,     ONLY: LSUPSAT
USE MODD_PARAMETERS,     ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_ICE_n,    ONLY: CSEDIM, LADJ_BEFORE, LADJ_AFTER, LRED, PARAM_ICEN
USE MODD_PARAM_LIMA,     ONLY: LADJ, LPTSPLIT, LSPRO, NMOD_CCN, NMOD_IFN, NMOD_IMM, NMOM_I
USE MODD_RAIN_ICE_DESCR_n, ONLY: XRTMIN, RAIN_ICE_DESCRN
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAMN
USE MODD_SALT,           ONLY: LSALT
USE MODD_TURB_n,         ONLY: TURBN
!
USE MODE_ll
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
use mode_sources_neg_correct, only: Sources_neg_correct
!
USE MODI_C2R2_ADJUST
USE MODI_FAST_TERMS
USE MODI_GET_HALO
USE MODI_ICE_ADJUST
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
USE MODI_RAIN_ICE_OLD
USE MODI_SHUMAN
USE MODI_SLOW_TERMS
USE MODI_AER2LIMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud paramerization
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
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE   ! Output file
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
real, dimension(:,:,:), allocatable :: ZEXN
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
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZICEFR
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZPRCFR
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZTM
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSIGQSAT2D
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZDUM
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
CALL FILL_DIMPHYEX(YLDIMPHYEX, SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3))
!
GWEST  = LWEST_ll()
GEAST  = LEAST_ll()
GSOUTH = LSOUTH_ll()
GNORTH = LNORTH_ll()
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

DEALLOCATE(ZSVT)
END IF

!UPG*PT
!
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
PTHS(:,:,:) = PTHS(:,:,:) / PRHODJ(:,:,:)
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
  DO JSV = ISVBEG, ISVEND
    PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) / PRHODJ(:,:,:)
  ENDDO
ENDIF
!
!  complete the lateral boundaries to avoid possible problems
!
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
IF(GWEST  .AND. HLBCX(1) /= 'CYCL')  PRT(:IIB-1,:,:,2:) = 0.0
IF(GEAST  .AND. HLBCX(2) /= 'CYCL')  PRT(IIE+1:,:,:,2:) = 0.0
IF(GSOUTH .AND. HLBCY(1) /= 'CYCL')  PRT(:,:IJB-1,:,2:) = 0.0
IF(GNORTH .AND. HLBCY(2) /= 'CYCL')  PRT(:,IJE+1:,:,2:) = 0.0
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
DO JI=1,JPHEXT
  PSVS(JI,     :,      :, ISVBEG:ISVEND) = PSVS(IIB, :,   :, ISVBEG:ISVEND)
  PSVS(IIE+JI, :,      :, ISVBEG:ISVEND) = PSVS(IIE, :,   :, ISVBEG:ISVEND)
  PSVS(:,      JI,     :, ISVBEG:ISVEND) = PSVS(:,   IJB, :, ISVBEG:ISVEND)
  PSVS(:,      IJE+JI, :, ISVBEG:ISVEND) = PSVS(:,   IJE, :, ISVBEG:ISVEND)
END DO
 !
!  complete the physical boundaries to avoid some computations
!
  IF(GWEST  .AND. HLBCX(1) /= 'CYCL') PSVT(:IIB-1, :,      :, ISVBEG:ISVEND) = 0.0
  IF(GEAST  .AND. HLBCX(2) /= 'CYCL') PSVT(IIE+1:, :,      :, ISVBEG:ISVEND) = 0.0
  IF(GSOUTH .AND. HLBCY(1) /= 'CYCL') PSVT(:,      :IJB-1, :, ISVBEG:ISVEND) = 0.0
  IF(GNORTH .AND. HLBCY(2) /= 'CYCL') PSVT(:,      IJE+1:, :, ISVBEG:ISVEND) = 0.0
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
!
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
! IF (      HCLOUD == 'KESS'                       &
!      .OR. HCLOUD == 'ICE3' .OR. HCLOUD == 'ICE4' &
!      .OR. HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' &
!      .OR. HCLOUD == 'KHKO' .OR. HCLOUD == 'LIMA' ) THEN
! !
!   DO JRR = 3,KRR
!     SELECT CASE (JRR)
!       CASE(3,5,6,7) ! rain, snow, graupel and hail
!
!         IF ( MIN_ll( PRS(:,:,:,JRR), IINFO_ll) < 0.0 ) THEN
! !
! ! compute the total water mass computation
! !
!           ZMASSTOT = MAX( 0. , SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
! !
! ! remove the negative values
! !
!           PRS(:,:,:,JRR) = MAX( 0., PRS(:,:,:,JRR) )
! !
! ! compute the new total mass
! !
!           ZMASSPOS = MAX(XMNH_TINY,SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
! !
! ! correct again in such a way to conserve the total mass
! !
!           ZRATIO = ZMASSTOT / ZMASSPOS
!           PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * ZRATIO
! !
!         END IF
!     END SELECT
!   END DO
! END IF
!
!*       3.2    Adjustement for liquid and solid cloud
!
! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
call Sources_neg_correct( hcloud, 'NEGA', krr, ptstep, ppabst, ptht, prt, pths, prs, psvs, prhodj )
!
!*       3.4    Limitations of Na and Nc to the CCN max number concentration
!
! Commented by O.Thouron 03/2013
!IF ((HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO') &
!     .AND.(XCONC_CCN > 0)) THEN
!  IF ((HACTCCN /= 'ABRK')) THEN
!  ZSVT(:,:,:,1) = MIN( ZSVT(:,:,:,1),XCONC_CCN )
!  ZSVT(:,:,:,2) = MIN( ZSVT(:,:,:,2),XCONC_CCN )
!  ZSVS(:,:,:,1) = MIN( ZSVS(:,:,:,1),XCONC_CCN )
!  ZSVS(:,:,:,2) = MIN( ZSVS(:,:,:,2),XCONC_CCN )
!  END IF
!END IF
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
    allocate( zexn( size( pzz, 1 ), size( pzz, 2 ), size( pzz, 3 ) ) )
    ZEXN(:,:,:)= (PPABST(:,:,:)/CST%XP00)**(CST%XRD/CST%XCPD)
!
!*       9.1    Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF( PZZ )
    IF(LRED .AND. LADJ_BEFORE) THEN
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,           &
                      PARAM_ICEN, TBUCONF, KRR,                                &
                      'ADJU',                                                  &
                      PTSTEP, ZSIGQSAT2D,                                      &
                      PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV,PMFCONV, PPABST, ZZZ,  &
                      ZEXN, PCF_MF, PRC_MF, PRI_MF,                            &
                      ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                            &
                      PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,        &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTH=PTHS*PTSTEP, PTHS=PTHS,                              &
                      OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR,    &
                      PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                      PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                      PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                      PRG=PRS(:,:,:,6)*PTSTEP,                                 &
                      TBUDGETS=TBUDGETS,KBUDGETS=SIZE(TBUDGETS),               &
                      PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                    &
                      PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF                     )
    ENDIF
    IF (LRED) THEN
      CALL RAIN_ICE (YLDIMPHYEX,CST, PARAM_ICEN, RAIN_ICE_PARAMN,        &
                    RAIN_ICE_DESCRN, TBUCONF,                            &
                    PTSTEP, KRR, ZEXN,                                   &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT,PCLDFR,&
                    PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,              &
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                    &
                    PRT(:,:,:,3), PRT(:,:,:,4),                          &
                    PRT(:,:,:,5), PRT(:,:,:,6),                          &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),      &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),            &
                    PINPRC,PINPRR, PEVAP3D,                    &
                    PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,              &
                    TBUDGETS,SIZE(TBUDGETS),           &
                    PSEA,PTOWN, PFPR=ZFPR                                )
    ELSE 
      CALL RAIN_ICE_OLD (YLDIMPHYEX, OSEDIC, CSEDIM, HSUBG_AUCV, OWARM, 1, IKU, 1,    &
                    KSPLITR, PTSTEP, KRR,                                 &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC,PINPRR, PINPRR3D, PEVAP3D,                     &
                    PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,                &
                    PSEA, PTOWN, PFPR=ZFPR)
    END IF

!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,           &
                       PARAM_ICEN, TBUCONF, KRR,                               &
                       'DEPI',                                                 &
                       PTSTEP, ZSIGQSAT2D,                                     &
                       PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,PPABST, ZZZ, &
                       ZEXN, PCF_MF, PRC_MF, PRI_MF,                           &
                       ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                           &
                       PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,       &
                       PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                   &
                       PTH=PTHS*PTSTEP, PTHS=PTHS,                             &
                       OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR, &
                       PRR=PRS(:,:,:,3)*PTSTEP,                                &
                       PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),             &
                       PRS=PRS(:,:,:,5)*PTSTEP,                                &
                       PRG=PRS(:,:,:,6)*PTSTEP,                                &
                       TBUDGETS=TBUDGETS,KBUDGETS=SIZE(TBUDGETS),              &
                       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                   &
                       PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF                    )
    END IF

    deallocate( zexn )
!
  CASE ('ICE4')
!
!*       10.    MIXED-PHASE MICROPHYSICAL SCHEME (WITH 4 ICE SPECIES)
!               -----------------------------------------------------
!
    allocate( zexn( size( pzz, 1 ), size( pzz, 2 ), size( pzz, 3 ) ) )
    ZEXN(:,:,:)= (PPABST(:,:,:)/CST%XP00)**(CST%XRD/CST%XCPD)
!
!*       10.1   Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF( PZZ )
    IF(LRED .AND. LADJ_BEFORE) THEN
      CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,           &
                       PARAM_ICEN, TBUCONF, KRR,                               &
                       'ADJU',                                                 &
                       PTSTEP, ZSIGQSAT2D,                                     &
                       PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV,PMFCONV, PPABST, ZZZ, &
                       ZEXN, PCF_MF, PRC_MF, PRI_MF,                           &
                       ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                           &
                       PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,       &
                       PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                   &
                       PTH=PTHS*PTSTEP, PTHS=PTHS,                             &
                       OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR, &
                       PRR=PRS(:,:,:,3)*PTSTEP,                                &
                       PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),             &
                       PRS=PRS(:,:,:,5)*PTSTEP,                                &
                       PRG=PRS(:,:,:,6)*PTSTEP,                                &
                       TBUDGETS=TBUDGETS,KBUDGETS=SIZE(TBUDGETS),              &
                       PRH=PRS(:,:,:,7)*PTSTEP,                                &
                       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                   &
                       PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF                    )
    ENDIF
    IF  (LRED) THEN
     CALL RAIN_ICE (YLDIMPHYEX,CST, PARAM_ICEN, RAIN_ICE_PARAMN,          &
                    RAIN_ICE_DESCRN, TBUCONF,                             &
                    PTSTEP, KRR, ZEXN,                                    &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC, PINPRR, PEVAP3D,                    &
                    PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,               &
                    TBUDGETS,SIZE(TBUDGETS),                     &            
                    PSEA, PTOWN,                                          &
                    PRT(:,:,:,7), PRS(:,:,:,7), PINPRH, PFPR=ZFPR         )
    ELSE
      CALL RAIN_ICE_OLD (YLDIMPHYEX, OSEDIC, CSEDIM, HSUBG_AUCV, OWARM, 1, IKU, 1,    &
                    KSPLITR, PTSTEP, KRR,                                 &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC,PINPRR, PINPRR3D, PEVAP3D,                     &
                    PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,                &
                    PSEA, PTOWN,                                          &
                    PRT(:,:,:,7),  PRS(:,:,:,7), PINPRH, PFPR=ZFPR)
    END IF


!
!*       10.2   Perform the saturation adjustment over cloud ice and cloud water
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
     CALL ICE_ADJUST (YLDIMPHYEX,CST, RAIN_ICE_PARAMN, NEBN, TURBN,          &
                     PARAM_ICEN, TBUCONF, KRR,                               &
                     'DEPI',                                                 &
                     PTSTEP, ZSIGQSAT2D,                                     &
                     PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,PPABST, ZZZ, &
                     ZEXN, PCF_MF, PRC_MF, PRI_MF,                           &
                     ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                           &
                     PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,       &
                     PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                   &
                     PTH=PTHS*PTSTEP, PTHS=PTHS,                             &
                     OCOMPUTE_SRC=SIZE(PSRCS, 3)/=0, PSRCS=PSRCS, PCLDFR=PCLDFR, &
                     PRR=PRS(:,:,:,3)*PTSTEP,                                &
                     PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),             &
                     PRS=PRS(:,:,:,5)*PTSTEP,                                &
                     PRG=PRS(:,:,:,6)*PTSTEP,                                &
                     TBUDGETS=TBUDGETS,KBUDGETS=SIZE(TBUDGETS),              &
                     PRH=PRS(:,:,:,7)*PTSTEP,                                &
                     PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF,                   &
                     PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF                    )
    END IF

    deallocate( zexn )
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
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF( PZZ )
     IF (LPTSPLIT) THEN
        CALL LIMA (YLDIMPHYEX,CST,TBUCONF,TBUDGETS,SIZE(TBUDGETS),         &
                   PTSTEP,                                                 &
                   PRHODREF, PEXNREF, ZDZZ,                                &
                   PRHODJ, PPABST,                                         &
                   NMOD_CCN, NMOD_IFN, NMOD_IMM,                           &
                   PDTHRAD, PTHT, PRT,                                     &
                   PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), PW_ACT,          &
                   PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),       &
                   PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, PINPRH, &
                   PEVAP3D, PCLDFR, PICEFR, PRAINFR, ZFPR                 )
     ELSE

        IF (OWARM) CALL LIMA_WARM(OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP, KMI,       &
                                  TPFILE, KRR, PZZ, PRHODJ,                         &
                                  PRHODREF, PEXNREF, PW_ACT, PPABST,                &
                                  PDTHRAD,                                          &
                                  PTHT, PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                  PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                                  PINPRC, PINPRR, PINDEP, PINPRR3D, PEVAP3D         )
!
        IF (NMOM_I.GE.1) CALL LIMA_COLD(CST, OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,    &
                                  KRR, PZZ, PRHODJ,                                  &
                                  PRHODREF, PEXNREF, PPABST, PW_ACT,                 &
                                  PTHT, PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),  &
                                  PTHS, PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),  &
                                  PINPRS, PINPRG, PINPRH                             )
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
   IF (LSPRO) THEN
    CALL LIMA_NOTADJUST (KMI, TPFILE, HRAD,                                      &
                         PTSTEP, PRHODJ, PPABSTT, PPABST, PRHODREF, PEXNREF, PZZ, &
                         PTHT,PRT, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),        &
                         PTHS,PRS, PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),        &
                         PCLDFR, PICEFR, PRAINFR, PSRCS                          )
   ELSE IF (LPTSPLIT) THEN
    CALL LIMA_ADJUST_SPLIT(YLDIMPHYEX,CST,TBUCONF,TBUDGETS,SIZE(TBUDGETS),           &
                     KRR, KMI, CCONDENS, CLAMBDA3,                                   &
                     OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,                          &
                     PRHODREF, PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, PPABSTT, ZZZ,&
                     PDTHRAD, PW_ACT,                                                &
                     PRT, PRS, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),                &
                     PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),                          &
                     PTHS, PSRCS, PCLDFR, PICEFR, PRC_MF, PRI_MF, PCF_MF             )
   ELSE
    CALL LIMA_ADJUST(KRR, KMI, TPFILE,                                &
                     OSUBG_COND, PTSTEP,                              &
                     PRHODREF, PRHODJ, PEXNREF, PPABST, PPABSTT,      &
                     PRT, PRS, PSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), &
                     PSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),           &
                     PTHS, PSRCS, PCLDFR, PICEFR, PRAINFR             )
   ENDIF
!
END SELECT
!
IF(HCLOUD=='ICE3' .OR. HCLOUD=='ICE4' ) THEN
! TODO: code a generic routine to update vertical lower and upper levels to 0, a
! specific value or to IKB or IKE and apply it to every output prognostic variable of physics
  PCIT(:,:,1)     = 0.
  PCIT(:,:,IKE+1) = 0.

  PINPRC3D=ZFPR(:,:,:,2) / CST%XRHOLW
  PINPRR3D=ZFPR(:,:,:,3) / CST%XRHOLW
  PINPRS3D=ZFPR(:,:,:,5) / CST%XRHOLW
  PINPRG3D=ZFPR(:,:,:,6) / CST%XRHOLW
  IF(KRR==7) PINPRH3D=ZFPR(:,:,:,7) / CST%XRHOLW
  WHERE (PRT(:,:,:,2) > 1.E-04 )
    PSPEEDC=ZFPR(:,:,:,2) / (PRT(:,:,:,2) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,3) > 1.E-04 )
    PSPEEDR=ZFPR(:,:,:,3) / (PRT(:,:,:,3) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,5) > 1.E-04 )
    PSPEEDS=ZFPR(:,:,:,5) / (PRT(:,:,:,5) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,6) > 1.E-04 )
    PSPEEDG=ZFPR(:,:,:,6) / (PRT(:,:,:,6) * PRHODREF(:,:,:))
  ENDWHERE
  IF(KRR==7) THEN
    WHERE (PRT(:,:,:,7) > 1.E-04 )
      PSPEEDH=ZFPR(:,:,:,7) / (PRT(:,:,:,7) * PRHODREF(:,:,:))
    ENDWHERE
  ENDIF
ENDIF

! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
call Sources_neg_correct( hcloud, 'NECON', krr, ptstep, ppabst, ptht, prt, pths, prs, psvs, prhodj )

!-------------------------------------------------------------------------------
!
!
!*      13.     SWITCH BACK TO THE PROGNOSTIC VARIABLES
!               ---------------------------------------
!
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

!-------------------------------------------------------------------------------
!
END SUBROUTINE RESOLVED_CLOUD
