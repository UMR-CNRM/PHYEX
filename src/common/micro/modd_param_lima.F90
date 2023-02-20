!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ######################
      MODULE MODD_PARAM_LIMA
!     ######################
!
!!****  *MODD_PARAM_LIMA* - declaration of the control parameters
!!                               for use in the LIMA scheme.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare the microphysical
!!    constants. This includes the descriptive parameters for the raindrop 
!!    and the parameters relevant of the dimensional distributions.
!!
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty  *Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe            14/03/2022  add CIBU and RDSF
!!
!-------------------------------------------------------------------------------
!
USE MODD_PARAMETERS, ONLY : JPLIMACCNMAX, JPLIMAIFNMAX
!
IMPLICIT NONE
!
LOGICAL, SAVE :: LLIMA_DIAG             ! Compute diagnostics for concentration /m3
!
LOGICAL, SAVE :: LPTSPLIT               ! activate time-splitting technique by S. Riette
LOGICAL, SAVE :: LFEEDBACKT             ! recompute tendencies if T changes sign
INTEGER, SAVE :: NMAXITER               ! maximum number of iterations
REAL,    SAVE :: XMRSTEP                ! maximum change in mixing ratio allowed before recomputing tedencies
REAL,    SAVE :: XTSTEP_TS              ! maximum time for the sub-time-step
!
!*       1.   COLD SCHEME
!             -----------
!
! 1.1 Cold scheme configuration
!
LOGICAL, SAVE :: LNUCL                  ! TRUE to enable ice nucleation
LOGICAL, SAVE :: LSEDI                  ! TRUE to enable pristine ice sedimentation
LOGICAL, SAVE :: LHHONI                 ! TRUE to enable freezing of haze particules
LOGICAL, SAVE :: LMEYERS                ! TRUE to use Meyers nucleation
LOGICAL, SAVE :: LCIBU                  ! TRUE to use collisional ice breakup
LOGICAL, SAVE :: LRDSF                  ! TRUE to use rain drop shattering by freezing
INTEGER, SAVE :: NMOM_I                 ! Number of moments for pristine ice
INTEGER, SAVE :: NMOM_S                 ! Number of moments for snow
INTEGER, SAVE :: NMOM_G                 ! Number of moments for graupel
INTEGER, SAVE :: NMOM_H                 ! Number of moments for hail
!
! 1.2 IFN initialisation
!
INTEGER, SAVE          :: NMOD_IFN               ! Number of IFN modes
REAL, DIMENSION(JPLIMAIFNMAX), SAVE :: XIFN_CONC ! Ref. concentration of IFN(#/L)
LOGICAL, SAVE          :: LIFN_HOM               ! True for z-homogeneous IFN concentrations
CHARACTER(LEN=8), SAVE :: CIFN_SPECIES           ! Internal mixing species definitions
CHARACTER(LEN=8), SAVE :: CINT_MIXING            ! Internal mixing type selection (pure DM1 ...)
INTEGER, SAVE          :: NMOD_IMM               ! Number of CCN modes acting by immersion
INTEGER, SAVE          :: NIND_SPECIE            ! CCN acting by immersion are considered pure
                                                 ! IFN of either DM = 1, BC = 2 or O = 3
INTEGER, DIMENSION(:), SAVE, ALLOCATABLE :: NIMM            ! Link between CCN and IMM modes
INTEGER, DIMENSION(:), SAVE, ALLOCATABLE :: NINDICE_CCN_IMM ! ??????????
INTEGER, SAVE                            :: NSPECIE         ! Internal mixing number of species
REAL, DIMENSION(:),    SAVE, ALLOCATABLE :: XMDIAM_IFN      ! Mean diameter of IFN modes
REAL, DIMENSION(:),    SAVE, ALLOCATABLE :: XSIGMA_IFN      ! Sigma of IFN modes
REAL, DIMENSION(:),    SAVE, ALLOCATABLE :: XRHO_IFN        ! Density of IFN modes 
REAL, DIMENSION(:,:),  SAVE, ALLOCATABLE :: XFRAC           ! Composition of each IFN mode
REAL, DIMENSION(:),    SAVE, ALLOCATABLE :: XFRAC_REF       ! AP compostion in Phillips 08
!
! 1.3 Ice characteristics
!
LOGICAL, SAVE :: LSNOW_T                     ! TRUE to enable snow param. after Wurtz 2021
LOGICAL, SAVE :: LMURAKAMI                   ! snow + liq -> graupel after Murakami (as in RAIN_ICE_RED)
CHARACTER(LEN=4), SAVE :: CPRISTINE_ICE_LIMA ! Pristine type PLAT, COLU or BURO
CHARACTER(LEN=4), SAVE :: CHEVRIMED_ICE_LIMA ! Heavily rimed type GRAU or HAIL
REAL,SAVE              :: XALPHAI,XNUI,    & ! Pristine ice   distribution parameters
                          XALPHAS,XNUS,    & ! Snow/aggregate distribution parameters
                          XALPHAG,XNUG       ! Graupel        distribution parameters
!
! 1.4 Phillips (2013) nucleation parameterization
!
INTEGER, SAVE          :: NPHILLIPS     ! =8 for Phillips08, =13 for Phillips13
!
REAL, DIMENSION(4), SAVE   :: XT0       ! Threshold of T in H_X for X={DM1,DM2,BC,O} [K]
REAL, DIMENSION(4), SAVE   :: XDT0      ! Range in T for transition of H_X near XT0 [K]
REAL, DIMENSION(4), SAVE   :: XDSI0     ! Range in Si for transition of H_X near XSI0
REAL,               SAVE   :: XSW0      ! Threshold of Sw in H_X 
REAL,               SAVE   :: XRHO_CFDC ! Air density at which CFDC data were reported [kg m**3]
REAL, DIMENSION(4), SAVE   :: XH        ! Fraction<<1 of aerosol for X={DM,BC,O}
REAL, DIMENSION(4), SAVE   :: XAREA1    ! Total surface of all aerosols in group X with
                            ! diameters between 0.1 and 1 µm, for X={DM1,DM2,BC,O} [m**2 kg**-1]
REAL,               SAVE   :: XGAMMA    ! Factor boosting IN concentration due to 
                                        ! bulk-liquid modes
!
REAL, DIMENSION(4), SAVE   :: XTX1      ! Threshold of T in Xi for X={DM1,DM2,BC,O} [K]
REAL, DIMENSION(4), SAVE   :: XTX2      ! Threshold of T in Xi for X={DM1,DM2,BC,O} [K]
!
REAL,DIMENSION(:), SAVE, ALLOCATABLE :: XABSCISS, XWEIGHT ! Gauss quadrature method 
INTEGER,           SAVE              :: NDIAM             ! Gauss quadrature accuracy 
!
! 1.5 Meyers (1992) nucleation parameterization
!
REAL,SAVE :: XFACTNUC_DEP,XFACTNUC_CON  ! Amplification factor for IN conc.
                                        !   DEP refers to DEPosition mode
                                        !   CON refers to CONtact    mode
!
! 1.6 Collisional Ice Break Up parameterization
!
REAL,SAVE :: XNDEBRIS_CIBU              ! Number of ice crystal debris produced
                                        ! by the break up of aggregate particles
!
!-------------------------------------------------------------------------------
!
!
!*       2.   WARM SCHEME
!             -----------
!
! 2.1 Warm scheme configuration
!
LOGICAL, SAVE :: LACTI         ! TRUE to enable CCN activation
LOGICAL, SAVE :: LSEDC         ! TRUE to enable the droplet sedimentation
LOGICAL, SAVE :: LACTIT        ! TRUE to enable the usage of dT/dt in CCN activation
LOGICAL, SAVE :: LBOUND        ! TRUE to enable the continuously replenishing
                               ! aerosol concentrations through the open
                               ! lateral boundaries -> boundaries.f90
LOGICAL, SAVE :: LDEPOC        ! Deposition of rc at 1st level above ground
LOGICAL, SAVE :: LACTTKE       ! TRUE to take into account TKE in W for activation
LOGICAL, SAVE :: LADJ          ! TRUE for adjustment procedure + Smax (false for diagnostic supersaturation)
LOGICAL, SAVE :: LSPRO         ! TRUE for prognostic supersaturation                     
LOGICAL, SAVE :: LKHKO         ! TRUE for Scu simulation (replicates the previous KHKO scheme)                     
LOGICAL, SAVE :: LKESSLERAC    ! TRUE for Kessler autoconversion (if NMOM_C=1)
!
INTEGER, SAVE :: NMOM_C        ! Number of moments for cloud droplets
INTEGER, SAVE :: NMOM_R        ! Number of moments for rain drops
!
! 2.2 CCN initialisation
!
INTEGER,         SAVE                     :: NMOD_CCN         ! Number of CCN modes
REAL, DIMENSION(JPLIMACCNMAX), SAVE       :: XCCN_CONC        ! CCN conc.  (#/cm3)
LOGICAL,         SAVE                     :: LCCN_HOM         ! True for z-homogeneous CCN concentrations
CHARACTER(LEN=8),SAVE                     :: CCCN_MODES       ! CCN modes characteristics (Jungfraujoch ...)
REAL, DIMENSION(:), SAVE, ALLOCATABLE     :: XR_MEAN_CCN,   & ! Mean radius of CCN modes
                                             XLOGSIG_CCN,   & ! Log of geometric dispersion of the CCN modes
                                             XRHO_CCN         ! Density of the CCN modes
REAL, DIMENSION(:), SAVE, ALLOCATABLE     :: XKHEN_MULTI,   & ! Parameters defining the CCN activation
                                             XMUHEN_MULTI,  & ! spectra for a multimodal aerosol distribution
                                             XBETAHEN_MULTI   ! 
REAL, DIMENSION(:,:,:) ,SAVE, ALLOCATABLE :: XCONC_CCN_TOT    ! Total aerosol number concentration
REAL, DIMENSION(:),     SAVE, ALLOCATABLE :: XLIMIT_FACTOR    ! compute CHEN ????????????
!
! 2.3 Water particles characteristics
!
REAL,SAVE     :: XALPHAR,XNUR,       & ! Raindrop      distribution parameters
                 XALPHAC,XNUC          ! Cloud droplet distribution parameters
!
! 2.4 CCN activation
!
CHARACTER(LEN=3),SAVE :: HPARAM_CCN = 'CPB'   ! Parameterization of the CCN activation
CHARACTER(LEN=3),SAVE :: HINI_CCN             ! Initialization type of CCN activation
CHARACTER(LEN=10),DIMENSION(JPLIMACCNMAX),SAVE :: HTYPE_CCN ! 'M' or 'C' CCN type
REAL,SAVE             :: XFSOLUB_CCN,       & ! Fractionnal solubility of the CCN
                         XACTEMP_CCN,       & ! Expected temperature of CCN activation
                         XAERDIFF, XAERHEIGHT ! For the vertical gradient of aerosol distribution
!
! Cloud droplet deposition
!
REAL, SAVE :: XVDEPOC
!
!-------------------------------------------------------------------------------
!
!
!*       3.   BELOW CLOUD SCAVENGING
!             ----------------------
!
LOGICAL, SAVE :: LSCAV           ! TRUE for aerosol scavenging by precipitations 
LOGICAL, SAVE :: LAERO_MASS      ! TRUE to compute the total aerosol mass scavenging rate 
!
INTEGER       :: NDIAMR = 20     ! Max Number of droplet for quadrature method  
INTEGER       :: NDIAMP = 20     ! Max Number of aerosol particle for quadrature method  
!
REAL, SAVE    :: XT0SCAV = 293.15  ! [K]
REAL, SAVE    :: XTREF = 273.15    ! [K]
REAL, SAVE    :: XNDO = 8.*1.0E6   ! [/m**4]
!
!-------------------------------------------------------------------------------
!
!
!*       4.   ATMOSPHERIC & OTHER PARAMETERS
!             ------------------------------
!
REAL, SAVE    :: XMUA0     = 1.711E-05  ![Pa.s] Air Viscosity at T=273.15K
REAL, SAVE    :: XT_SUTH_A = 110.4      ![K] Sutherland Temperature for Air
REAL, SAVE    :: XMFPA0    = 6.6E-08    ![m] Mean Free Path of Air under standard conditions
!
REAL, SAVE    :: XVISCW = 1.0E-3        ![Pa.s] water viscosity at 20°C
! Correction
!REAL, SAVE    :: XRHO00 = 1.292        !rho on the floor    [Kg/m**3]
REAL, SAVE    :: XRHO00 = 1.2041        !rho at P=1013.25 and T=20°C
!
REAL,SAVE :: XCEXVT                     ! air density fall speed correction
!
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XRTMIN ! Min values of the mixing ratios
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XCTMIN ! Min values of the drop concentrations
!
!
! Sedimentation variables
!
INTEGER,DIMENSION(7),SAVE :: NSPLITSED
REAL,DIMENSION(7),SAVE :: XLB
REAL,DIMENSION(7),SAVE :: XLBEX
REAL,DIMENSION(7),SAVE :: XD
REAL,DIMENSION(7),SAVE :: XFSEDR
REAL,DIMENSION(7),SAVE :: XFSEDC
!
END MODULE MODD_PARAM_LIMA
