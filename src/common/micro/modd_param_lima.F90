!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ######################
      MODULE MODD_PARAM_LIMA
!     ######################
!> @file
!!      *MODD_PARAM_LIMA* - declaration of the control parameters
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
!!    IMPLICIT ARGUMENTS
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
TYPE PARAM_LIMA_t
LOGICAL :: LLIMA_DIAG             ! Compute diagnostics for concentration /m3
!
LOGICAL :: LPTSPLIT               ! activate time-splitting technique by S. Riette
LOGICAL :: LFEEDBACKT             ! recompute tendencies if T changes sign
INTEGER :: NMAXITER               ! maximum number of iterations
REAL    :: XMRSTEP                ! maximum change in mixing ratio allowed before recomputing tedencies
REAL    :: XTSTEP_TS              ! maximum time for the sub-time-step
!
!*       1.   COLD SCHEME
!             -----------
!
! 1.1 Cold scheme configuration
!
LOGICAL :: LNUCL                  ! TRUE to enable ice nucleation
LOGICAL :: LSEDI                  ! TRUE to enable pristine ice sedimentation
LOGICAL :: LHHONI                 ! TRUE to enable freezing of haze particules
LOGICAL :: LMEYERS                ! TRUE to use Meyers nucleation
LOGICAL :: LCIBU                  ! TRUE to use collisional ice breakup
LOGICAL :: LRDSF                  ! TRUE to use rain drop shattering by freezing
INTEGER :: NMOM_I                 ! Number of moments for pristine ice
INTEGER :: NMOM_S                 ! Number of moments for snow
INTEGER :: NMOM_G                 ! Number of moments for graupel
INTEGER :: NMOM_H                 ! Number of moments for hail
!
! 1.2 IFN initialisation
!
INTEGER          :: NMOD_IFN               ! Number of IFN modes
REAL, DIMENSION(JPLIMAIFNMAX) :: XIFN_CONC ! Ref. concentration of IFN(#/L)
LOGICAL          :: LIFN_HOM               ! True for z-homogeneous IFN concentrations
CHARACTER(LEN=8) :: CIFN_SPECIES           ! Internal mixing species definitions
CHARACTER(LEN=8) :: CINT_MIXING            ! Internal mixing type selection (pure DM1 ...)
INTEGER          :: NMOD_IMM               ! Number of CCN modes acting by immersion
INTEGER          :: NIND_SPECIE            ! CCN acting by immersion are considered pure
                                           ! IFN of either DM = 1, BC = 2 or O = 3
INTEGER, DIMENSION(:), ALLOCATABLE :: NIMM            ! Link between CCN and IMM modes
INTEGER, DIMENSION(:), ALLOCATABLE :: NINDICE_CCN_IMM ! ??????????
INTEGER                            :: NSPECIE         ! Internal mixing number of species
REAL, DIMENSION(:),    ALLOCATABLE :: XMDIAM_IFN      ! Mean diameter of IFN modes
REAL, DIMENSION(:),    ALLOCATABLE :: XSIGMA_IFN      ! Sigma of IFN modes
REAL, DIMENSION(:),    ALLOCATABLE :: XRHO_IFN        ! Density of IFN modes 
REAL, DIMENSION(:,:),  ALLOCATABLE :: XFRAC           ! Composition of each IFN mode
REAL, DIMENSION(:),    ALLOCATABLE :: XFRAC_REF       ! AP compostion in Phillips 08
!
! 1.3 Ice characteristics
!
LOGICAL :: LSNOW_T                     ! TRUE to enable snow param. after Wurtz 2021
LOGICAL :: LMURAKAMI                   ! snow + liq -> graupel after Murakami (as in RAIN_ICE_RED)
CHARACTER(LEN=4) :: CPRISTINE_ICE_LIMA ! Pristine type PLAT, COLU or BURO
CHARACTER(LEN=4) :: CHEVRIMED_ICE_LIMA ! Heavily rimed type GRAU or HAIL
REAL                   :: XALPHAI,XNUI,    & ! Pristine ice   distribution parameters
                          XALPHAS,XNUS,    & ! Snow/aggregate distribution parameters
                          XALPHAG,XNUG       ! Graupel        distribution parameters
!
! 1.4 Phillips (2013) nucleation parameterization
!
INTEGER              :: NPHILLIPS     ! =8 for Phillips08, =13 for Phillips13
!
REAL, DIMENSION(4)   :: XT0       ! Threshold of T in H_X for X={DM1,DM2,BC,O} [K]
REAL, DIMENSION(4)   :: XDT0      ! Range in T for transition of H_X near XT0 [K]
REAL, DIMENSION(4)   :: XDSI0     ! Range in Si for transition of H_X near XSI0
REAL                 :: XSW0      ! Threshold of Sw in H_X 
REAL                 :: XRHO_CFDC ! Air density at which CFDC data were reported [kg m**3]
REAL, DIMENSION(4)   :: XH        ! Fraction<<1 of aerosol for X={DM,BC,O}
REAL, DIMENSION(4)   :: XAREA1    ! Total surface of all aerosols in group X with
                                  ! diameters between 0.1 and 1 µm, for X={DM1,DM2,BC,O} [m**2 kg**-1]
REAL                 :: XGAMMA    ! Factor boosting IN concentration due to 
                                        ! bulk-liquid modes
!
REAL, DIMENSION(4)   :: XTX1      ! Threshold of T in Xi for X={DM1,DM2,BC,O} [K]
REAL, DIMENSION(4)   :: XTX2      ! Threshold of T in Xi for X={DM1,DM2,BC,O} [K]
!
REAL,DIMENSION(:), ALLOCATABLE :: XABSCISS, XWEIGHT ! Gauss quadrature method 
INTEGER                        :: NDIAM             ! Gauss quadrature accuracy 
!
! 1.5 Meyers (1992) nucleation parameterization
!
REAL      :: XFACTNUC_DEP,XFACTNUC_CON  ! Amplification factor for IN conc.
                                        !   DEP refers to DEPosition mode
                                        !   CON refers to CONtact    mode
!
! 1.6 Collisional Ice Break Up parameterization
!
REAL      :: XNDEBRIS_CIBU              ! Number of ice crystal debris produced
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
LOGICAL :: LACTI         ! TRUE to enable CCN activation
LOGICAL :: LSEDC         ! TRUE to enable the droplet sedimentation
LOGICAL :: LACTIT        ! TRUE to enable the usage of dT/dt in CCN activation
LOGICAL :: LDEPOC        ! Deposition of rc at 1st level above ground
LOGICAL :: LACTTKE       ! TRUE to take into account TKE in W for activation
LOGICAL :: LADJ          ! TRUE for adjustment procedure + Smax (false for diagnostic supersaturation)
LOGICAL :: LSPRO         ! TRUE for prognostic supersaturation                     
LOGICAL :: LKHKO         ! TRUE for Scu simulation (replicates the previous KHKO scheme)                     
LOGICAL :: LKESSLERAC    ! TRUE for Kessler autoconversion (if NMOM_C=1)
!
INTEGER :: NMOM_C        ! Number of moments for cloud droplets
INTEGER :: NMOM_R        ! Number of moments for rain drops
!
! 2.2 CCN initialisation
!
INTEGER                             :: NMOD_CCN         ! Number of CCN modes
REAL, DIMENSION(JPLIMACCNMAX)       :: XCCN_CONC        ! CCN conc.  (#/cm3)
LOGICAL                             :: LCCN_HOM         ! True for z-homogeneous CCN concentrations
CHARACTER(LEN=8)                    :: CCCN_MODES       ! CCN modes characteristics (Jungfraujoch ...)
REAL, DIMENSION(:), ALLOCATABLE     :: XR_MEAN_CCN,   & ! Mean radius of CCN modes
                                       XLOGSIG_CCN,   & ! Log of geometric dispersion of the CCN modes
                                       XRHO_CCN         ! Density of the CCN modes
REAL, DIMENSION(:), ALLOCATABLE     :: XKHEN_MULTI,   & ! Parameters defining the CCN activation
                                       XMUHEN_MULTI,  & ! spectra for a multimodal aerosol distribution
                                       XBETAHEN_MULTI   ! 
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCONC_CCN_TOT    ! Total aerosol number concentration
REAL, DIMENSION(:),     ALLOCATABLE :: XLIMIT_FACTOR    ! compute CHEN ????????????
!
! 2.3 Water particles characteristics
!
REAL          :: XALPHAR,XNUR,       & ! Raindrop      distribution parameters
                 XALPHAC,XNUC          ! Cloud droplet distribution parameters
!
! 2.4 CCN activation
!
CHARACTER(LEN=3)      :: HPARAM_CCN = 'CPB'   ! Parameterization of the CCN activation
CHARACTER(LEN=3)      :: HINI_CCN             ! Initialization type of CCN activation
CHARACTER(LEN=10), DIMENSION(JPLIMACCNMAX) :: HTYPE_CCN ! 'M' or 'C' CCN type
REAL                  :: XFSOLUB_CCN,       & ! Fractionnal solubility of the CCN
                         XACTEMP_CCN,       & ! Expected temperature of CCN activation
                         XAERDIFF, XAERHEIGHT ! For the vertical gradient of aerosol distribution
!
! Cloud droplet deposition
!
REAL :: XVDEPOC
!
!-------------------------------------------------------------------------------
!
!
!*       3.   BELOW CLOUD SCAVENGING
!             ----------------------
!
LOGICAL :: LSCAV           ! TRUE for aerosol scavenging by precipitations 
LOGICAL :: LAERO_MASS      ! TRUE to compute the total aerosol mass scavenging rate 
!
INTEGER :: NDIAMR = 20     ! Max Number of droplet for quadrature method  
INTEGER :: NDIAMP = 20     ! Max Number of aerosol particle for quadrature method  
!
REAL    :: XT0SCAV = 293.15  ! [K]
REAL    :: XTREF = 273.15    ! [K]
REAL    :: XNDO = 8.*1.0E6   ! [/m**4]
!
!-------------------------------------------------------------------------------
!
!
!*       4.   ATMOSPHERIC & OTHER PARAMETERS
!             ------------------------------
!
REAL    :: XMUA0     = 1.711E-05  ![Pa.s] Air Viscosity at T=273.15K
REAL    :: XT_SUTH_A = 110.4      ![K] Sutherland Temperature for Air
REAL    :: XMFPA0    = 6.6E-08    ![m] Mean Free Path of Air under standard conditions
!
REAL    :: XVISCW = 1.0E-3        ![Pa.s] water viscosity at 20°C
! Correction
!REAL    :: XRHO00 = 1.292        !rho on the floor    [Kg/m**3]
REAL    :: XRHO00 = 1.2041        !rho at P=1013.25 and T=20°C
!
REAL    :: XCEXVT                     ! air density fall speed correction
!
REAL, DIMENSION(:), ALLOCATABLE :: XRTMIN ! Min values of the mixing ratios
REAL, DIMENSION(:), ALLOCATABLE :: XCTMIN ! Min values of the drop concentrations
!
!
! Sedimentation variables
!
INTEGER,DIMENSION(7)   :: NSPLITSED
REAL,DIMENSION(7)      :: XLB
REAL,DIMENSION(7)      :: XLBEX
REAL,DIMENSION(7)      :: XD
REAL,DIMENSION(7)      :: XFSEDR
REAL,DIMENSION(7)      :: XFSEDC
END TYPE PARAM_LIMA_t
!
TYPE(PARAM_LIMA_t), TARGET, SAVE :: PARAM_LIMA
!
LOGICAL, POINTER :: LLIMA_DIAG => NULL(), &
                    LPTSPLIT => NULL(), &
                    LFEEDBACKT => NULL(), &
                    LNUCL => NULL(), &
                    LSEDI => NULL(), &
                    LHHONI => NULL(), &
                    LMEYERS => NULL(), &
                    LCIBU => NULL(), &
                    LRDSF => NULL(), &
                    LIFN_HOM => NULL(), &
                    LSNOW_T => NULL(), &
                    LMURAKAMI => NULL(), &
                    LACTI => NULL(), &
                    LSEDC => NULL(), &
                    LACTIT => NULL(), &
                    LDEPOC => NULL(), &
                    LACTTKE => NULL(), &
                    LADJ => NULL(), &
                    LSPRO => NULL(), &
                    LKHKO => NULL(), &
                    LKESSLERAC => NULL(), &
                    LCCN_HOM => NULL(), &
                    LSCAV => NULL(), &
                    LAERO_MASS => NULL()

INTEGER, POINTER :: NMAXITER => NULL(), &
                    NMOM_I => NULL(), &
                    NMOM_S => NULL(), &
                    NMOM_G => NULL(), &
                    NMOM_H => NULL(), &
                    NMOD_IFN => NULL(), &
                    NMOD_IMM => NULL(), &
                    NIND_SPECIE => NULL(), &
                    NSPECIE => NULL(), &
                    NPHILLIPS => NULL(), &
                    NDIAM => NULL(), &
                    NMOM_C => NULL(), &
                    NMOM_R => NULL(), &
                    NMOD_CCN => NULL(), &
                    NDIAMR => NULL(), &
                    NDIAMP => NULL()

REAL, POINTER :: XMRSTEP => NULL(), &
                 XTSTEP_TS => NULL(), &
                 XALPHAI => NULL(), &
                 XNUI => NULL(), &
                 XALPHAS => NULL(), &
                 XNUS => NULL(), &
                 XALPHAG => NULL(), &
                 XNUG => NULL(), &
                 XSW0 => NULL(), &
                 XRHO_CFDC => NULL(), &
                 XGAMMA => NULL(), &
                 XFACTNUC_DEP => NULL(), &
                 XFACTNUC_CON => NULL(), &
                 XNDEBRIS_CIBU => NULL(), &
                 XALPHAR => NULL(), &
                 XNUR => NULL(), &
                 XALPHAC => NULL(), &
                 XNUC => NULL(), &
                 XFSOLUB_CCN => NULL(), &
                 XACTEMP_CCN => NULL(), &
                 XAERDIFF => NULL(), &
                 XAERHEIGHT => NULL(), &
                 XVDEPOC => NULL(), &
                 XT0SCAV => NULL(), &
                 XTREF => NULL(), &
                 XNDO => NULL(), &
                 XMUA0 => NULL(), &
                 XT_SUTH_A => NULL(), &
                 XMFPA0 => NULL(), &
                 XVISCW => NULL(), &
                 XRHO00 => NULL(), &
                 XCEXVT => NULL()

REAL, DIMENSION(:), POINTER :: XIFN_CONC => NULL(), &
                               XMDIAM_IFN => NULL(), &
                               XSIGMA_IFN => NULL(), &
                               XRHO_IFN => NULL(), &
                               XFRAC_REF => NULL(), &
                               XT0 => NULL(), &
                               XDT0 => NULL(), &
                               XDSI0 => NULL(), &
                               XH => NULL(), &
                               XAREA1 => NULL(), &
                               XTX1 => NULL(), &
                               XTX2 => NULL(), &
                               XABSCISS => NULL(), &
                               XWEIGHT => NULL(), &
                               XCCN_CONC => NULL(), &
                               XR_MEAN_CCN => NULL(), &
                               XLOGSIG_CCN => NULL(), &
                               XRHO_CCN => NULL(), &
                               XKHEN_MULTI => NULL(), &
                               XMUHEN_MULTI => NULL(), &
                               XBETAHEN_MULTI => NULL(), &
                               XLIMIT_FACTOR => NULL(), &
                               XRTMIN => NULL(), &
                               XCTMIN => NULL(), &
                               XLB => NULL(), &
                               XLBEX => NULL(), &
                               XD => NULL(), &
                               XFSEDR => NULL(), &
                               XFSEDC => NULL()

REAL, DIMENSION(:,:),  POINTER :: XFRAC => NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCONC_CCN_TOT => NULL()

INTEGER, DIMENSION(:), POINTER :: NIMM => NULL(), &
                                  NINDICE_CCN_IMM => NULL(), &
                                  NSPLITSED => NULL()

CHARACTER(LEN=8), POINTER :: CIFN_SPECIES => NULL()
CHARACTER(LEN=8), POINTER :: CINT_MIXING => NULL()
CHARACTER(LEN=4), POINTER :: CPRISTINE_ICE_LIMA => NULL()
CHARACTER(LEN=4), POINTER :: CHEVRIMED_ICE_LIMA => NULL()
CHARACTER(LEN=8), POINTER :: CCCN_MODES => NULL()
CHARACTER(LEN=3), POINTER :: HPARAM_CCN => NULL()
CHARACTER(LEN=3), POINTER :: HINI_CCN => NULL()
CHARACTER(LEN=10), DIMENSION(:), POINTER :: HTYPE_CCN

NAMELIST/NAM_PARAM_LIMA/LNUCL, LSEDI, LHHONI, LMEYERS,                     &         
                        NMOM_I, NMOM_S, NMOM_G, NMOM_H,                    & 
                        NMOD_IFN, XIFN_CONC, LIFN_HOM,                     &
                        CIFN_SPECIES, CINT_MIXING, NMOD_IMM, NIND_SPECIE,  &
                        LSNOW_T, CPRISTINE_ICE_LIMA, CHEVRIMED_ICE_LIMA,   &                                   
                        !XALPHAI, XNUI, XALPHAS, XNUS, XALPHAG, XNUG,       &    
                        XFACTNUC_DEP, XFACTNUC_CON, NPHILLIPS,             &    
                        LCIBU, XNDEBRIS_CIBU, LRDSF, LMURAKAMI,            &                                         
                        LACTI, LSEDC, LACTIT, LSPRO,                       &                                         
                        LADJ, LKHKO, LKESSLERAC, NMOM_C, NMOM_R,           &                                         
                        NMOD_CCN, XCCN_CONC,                               &                                         
                        LCCN_HOM, CCCN_MODES, HINI_CCN, HTYPE_CCN,         &                                         
                        XALPHAC, XNUC, XALPHAR, XNUR,                      &                                         
                        XFSOLUB_CCN, XACTEMP_CCN, XAERDIFF, XAERHEIGHT,    &                                         
                        LSCAV, LAERO_MASS, LDEPOC, XVDEPOC, LACTTKE,       &                                         
                        LPTSPLIT, LFEEDBACKT, NMAXITER, XMRSTEP, XTSTEP_TS

CONTAINS
SUBROUTINE PARAM_LIMA_ASSOCIATE()
IMPLICIT NONE

IF(.NOT. ASSOCIATED(LLIMA_DIAG)) THEN
  LLIMA_DIAG         => PARAM_LIMA%LLIMA_DIAG          
  LPTSPLIT           => PARAM_LIMA%LPTSPLIT
  LFEEDBACKT         => PARAM_LIMA%LFEEDBACKT
  LNUCL              => PARAM_LIMA%LNUCL
  LSEDI              => PARAM_LIMA%LSEDI
  LHHONI             => PARAM_LIMA%LHHONI
  LMEYERS            => PARAM_LIMA%LMEYERS
  LCIBU              => PARAM_LIMA%LCIBU
  LRDSF              => PARAM_LIMA%LRDSF
  LIFN_HOM           => PARAM_LIMA%LIFN_HOM
  LSNOW_T            => PARAM_LIMA%LSNOW_T
  LMURAKAMI          => PARAM_LIMA%LMURAKAMI
  LACTI              => PARAM_LIMA%LACTI
  LSEDC              => PARAM_LIMA%LSEDC
  LACTIT             => PARAM_LIMA%LACTIT
  LDEPOC             => PARAM_LIMA%LDEPOC
  LACTTKE            => PARAM_LIMA%LACTTKE
  LADJ               => PARAM_LIMA%LADJ
  LSPRO              => PARAM_LIMA%LSPRO
  LKHKO              => PARAM_LIMA%LKHKO
  LKESSLERAC         => PARAM_LIMA%LKESSLERAC
  LCCN_HOM           => PARAM_LIMA%LCCN_HOM
  LSCAV              => PARAM_LIMA%LSCAV
  LAERO_MASS         => PARAM_LIMA%LAERO_MASS

  NMAXITER           => PARAM_LIMA%NMAXITER
  NMOM_I             => PARAM_LIMA%NMOM_I
  NMOM_S             => PARAM_LIMA%NMOM_S
  NMOM_G             => PARAM_LIMA%NMOM_G
  NMOM_H             => PARAM_LIMA%NMOM_H
  NMOD_IFN           => PARAM_LIMA%NMOD_IFN
  NMOD_IMM           => PARAM_LIMA%NMOD_IMM
  NIND_SPECIE        => PARAM_LIMA%NIND_SPECIE
  NSPECIE            => PARAM_LIMA%NSPECIE
  NPHILLIPS          => PARAM_LIMA%NPHILLIPS
  NDIAM              => PARAM_LIMA%NDIAM
  NMOM_C             => PARAM_LIMA%NMOM_C
  NMOM_R             => PARAM_LIMA%NMOM_R
  NMOD_CCN           => PARAM_LIMA%NMOD_CCN
  NDIAMR             => PARAM_LIMA%NDIAMR
  NDIAMP             => PARAM_LIMA%NDIAMP

  XMRSTEP            => PARAM_LIMA%XMRSTEP
  XTSTEP_TS          => PARAM_LIMA%XTSTEP_TS
  XALPHAI            => PARAM_LIMA%XALPHAI
  XNUI               => PARAM_LIMA%XNUI
  XALPHAS            => PARAM_LIMA%XALPHAS
  XNUS               => PARAM_LIMA%XNUS
  XALPHAG            => PARAM_LIMA%XALPHAG
  XNUG               => PARAM_LIMA%XNUG
  XSW0               => PARAM_LIMA%XSW0
  XRHO_CFDC          => PARAM_LIMA%XRHO_CFDC
  XGAMMA             => PARAM_LIMA%XGAMMA
  XFACTNUC_DEP       => PARAM_LIMA%XFACTNUC_DEP
  XFACTNUC_CON       => PARAM_LIMA%XFACTNUC_CON
  XNDEBRIS_CIBU      => PARAM_LIMA%XNDEBRIS_CIBU
  XALPHAR            => PARAM_LIMA%XALPHAR
  XNUR               => PARAM_LIMA%XNUR
  XALPHAC            => PARAM_LIMA%XALPHAC
  XNUC               => PARAM_LIMA%XNUC
  XFSOLUB_CCN        => PARAM_LIMA%XFSOLUB_CCN
  XACTEMP_CCN        => PARAM_LIMA%XACTEMP_CCN
  XAERDIFF           => PARAM_LIMA%XAERDIFF
  XAERHEIGHT         => PARAM_LIMA%XAERHEIGHT
  XVDEPOC            => PARAM_LIMA%XVDEPOC
  XT0SCAV            => PARAM_LIMA%XT0SCAV
  XTREF              => PARAM_LIMA%XTREF
  XNDO               => PARAM_LIMA%XNDO
  XMUA0              => PARAM_LIMA%XMUA0
  XT_SUTH_A          => PARAM_LIMA%XT_SUTH_A
  XMFPA0             => PARAM_LIMA%XMFPA0
  XVISCW             => PARAM_LIMA%XVISCW
  XRHO00             => PARAM_LIMA%XRHO00
  XCEXVT             => PARAM_LIMA%XCEXVT

  XIFN_CONC          => PARAM_LIMA%XIFN_CONC
  XT0                => PARAM_LIMA%XT0
  XDT0               => PARAM_LIMA%XDT0
  XDSI0              => PARAM_LIMA%XDSI0
  XH                 => PARAM_LIMA%XH
  XAREA1             => PARAM_LIMA%XAREA1
  XTX1               => PARAM_LIMA%XTX1
  XTX2               => PARAM_LIMA%XTX2
  XCCN_CONC          => PARAM_LIMA%XCCN_CONC
  XLB                => PARAM_LIMA%XLB
  XLBEX              => PARAM_LIMA%XLBEX
  XD                 => PARAM_LIMA%XD
  XFSEDR             => PARAM_LIMA%XFSEDR
  XFSEDC             => PARAM_LIMA%XFSEDC

  NSPLITSED          => PARAM_LIMA%NSPLITSED

  CIFN_SPECIES       => PARAM_LIMA%CIFN_SPECIES
  CINT_MIXING        => PARAM_LIMA%CINT_MIXING
  CPRISTINE_ICE_LIMA => PARAM_LIMA%CPRISTINE_ICE_LIMA
  CHEVRIMED_ICE_LIMA => PARAM_LIMA%CHEVRIMED_ICE_LIMA
  CCCN_MODES         => PARAM_LIMA%CCCN_MODES
  HPARAM_CCN         => PARAM_LIMA%HPARAM_CCN
  HINI_CCN           => PARAM_LIMA%HINI_CCN
  HTYPE_CCN          => PARAM_LIMA%HTYPE_CCN
ENDIF
END SUBROUTINE PARAM_LIMA_ASSOCIATE
!
SUBROUTINE PARAM_LIMA_DEALLOCATE(HNAME)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  SELECT CASE(TRIM(HNAME))
    CASE('NINDICE_CCN_IMM')
      DEALLOCATE(PARAM_LIMA%NINDICE_CCN_IMM)
      NINDICE_CCN_IMM => NULL()
  END SELECT
END SUBROUTINE PARAM_LIMA_DEALLOCATE
!
SUBROUTINE PARAM_LIMA_ALLOCATE(HNAME, KDIM1, KDIM2, KDIM3)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: HNAME
  INTEGER, INTENT(IN)           :: KDIM1
  INTEGER, OPTIONAL, INTENT(IN) :: KDIM2
  INTEGER, OPTIONAL, INTENT(IN) :: KDIM3

  SELECT CASE(TRIM(HNAME))
    !1D arrays
    CASE('NIMM')
      ALLOCATE(PARAM_LIMA%NIMM(KDIM1))
      NIMM => PARAM_LIMA%NIMM
    CASE('NINDICE_CCN_IMM')
      ALLOCATE(PARAM_LIMA%NINDICE_CCN_IMM(KDIM1))
      NINDICE_CCN_IMM => PARAM_LIMA%NINDICE_CCN_IMM
    CASE('XMDIAM_IFN')
      ALLOCATE(PARAM_LIMA%XMDIAM_IFN(KDIM1))
      XMDIAM_IFN => PARAM_LIMA%XMDIAM_IFN
    CASE('XSIGMA_IFN')
      ALLOCATE(PARAM_LIMA%XSIGMA_IFN(KDIM1))
      XSIGMA_IFN => PARAM_LIMA%XSIGMA_IFN
    CASE('XRHO_IFN')
      ALLOCATE(PARAM_LIMA%XRHO_IFN(KDIM1))
      XRHO_IFN => PARAM_LIMA%XRHO_IFN
    CASE('XFRAC_REF')
      ALLOCATE(PARAM_LIMA%XFRAC_REF(KDIM1))
      XFRAC_REF => PARAM_LIMA%XFRAC_REF
    CASE('XABSCISS')
      ALLOCATE(PARAM_LIMA%XABSCISS(KDIM1))
      XABSCISS => PARAM_LIMA%XABSCISS
    CASE('XWEIGHT')
      ALLOCATE(PARAM_LIMA%XWEIGHT(KDIM1))
      XWEIGHT => PARAM_LIMA%XWEIGHT
    CASE('XR_MEAN_CCN')
      ALLOCATE(PARAM_LIMA%XR_MEAN_CCN(KDIM1))
      XR_MEAN_CCN => PARAM_LIMA%XR_MEAN_CCN
    CASE('XLOGSIG_CCN')
      ALLOCATE(PARAM_LIMA%XLOGSIG_CCN(KDIM1))
      XLOGSIG_CCN => PARAM_LIMA%XLOGSIG_CCN
    CASE('XRHO_CCN')
      ALLOCATE(PARAM_LIMA%XRHO_CCN(KDIM1))
      XRHO_CCN => PARAM_LIMA%XRHO_CCN
    CASE('XKHEN_MULTI')
      ALLOCATE(PARAM_LIMA%XKHEN_MULTI(KDIM1))
      XKHEN_MULTI => PARAM_LIMA%XKHEN_MULTI
    CASE('XMUHEN_MULTI')
      ALLOCATE(PARAM_LIMA%XMUHEN_MULTI(KDIM1))
      XMUHEN_MULTI => PARAM_LIMA%XMUHEN_MULTI
    CASE('XBETAHEN_MULTI')
      ALLOCATE(PARAM_LIMA%XBETAHEN_MULTI(KDIM1))
      XBETAHEN_MULTI => PARAM_LIMA%XBETAHEN_MULTI
    CASE('XLIMIT_FACTOR')
      ALLOCATE(PARAM_LIMA%XLIMIT_FACTOR(KDIM1))
      XLIMIT_FACTOR => PARAM_LIMA%XLIMIT_FACTOR
    CASE('XRTMIN')
      ALLOCATE(PARAM_LIMA%XRTMIN(KDIM1))
      XRTMIN => PARAM_LIMA%XRTMIN
    CASE('XCTMIN')
      ALLOCATE(PARAM_LIMA%XCTMIN(KDIM1))
      XCTMIN => PARAM_LIMA%XCTMIN
    !
    !2D arrays
    CASE('XFRAC')
      ALLOCATE(PARAM_LIMA%XFRAC(KDIM1, KDIM2))
      XFRAC => PARAM_LIMA%XFRAC
    !
    !3D arrays
!    CASE('XCONC_CCN_TOT')
!      ALLOCATE(PARAM_LIMA%XCONC_CCN_TOT(KDIM1, KDIM2))
!      XCONC_CCN_TOT => PARAM_LIMA%XCONC_CCN_TOT
  END SELECT
END SUBROUTINE PARAM_LIMA_ALLOCATE
!
SUBROUTINE PARAM_LIMA_INIT(HPROGRAM, TFILENAM, LDNEEDNAM, KLUOUT, &
                          &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
!!*** *PARAM_ICEN_INIT* - Code needed to initialize the MODD_PARAM_LIMA module
!!
!!*   PURPOSE
!!    -------
!!    Sets the default values, reads the namelist, performs the checks and prints
!!
!!*   METHOD
!!    ------
!!    0. Declarations
!!       1. Declaration of arguments
!!       2. Declaration of local variables
!!    1. Default values
!!    2. Namelist
!!    3. Checks
!!    4. Prints
!!
!!    AUTHOR
!!    ------
!!    S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Apr 2023
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!       ---------------
!
USE MODE_POSNAM_PHY, ONLY: POSNAM_PHY
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_CHECK_NAM_VAL, ONLY: CHECK_NAM_VAL_CHAR
USE MODD_IO,  ONLY: TFILEDATA
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Name of the calling program
TYPE(TFILEDATA),   INTENT(IN) :: TFILENAM     !< Namelist file
LOGICAL,           INTENT(IN) :: LDNEEDNAM    !< True to abort if namelist is absent
INTEGER,           INTENT(IN) :: KLUOUT       !< Logical unit for outputs
LOGICAL, OPTIONAL, INTENT(IN) :: LDDEFAULTVAL !< Must we initialize variables with default values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDREADNAM    !< Must we read the namelist (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHECK      !< Must we perform some checks on values (defaults to .TRUE.)
INTEGER, OPTIONAL, INTENT(IN) :: KPRINT       !< Print level (defaults to 0): 0 for no print, 1 to safely print namelist,
                                              !! 2 to print informative messages
!
!* 0.2 Declaration of local variables
!      ------------------------------
!
LOGICAL :: LLDEFAULTVAL, LLREADNAM, LLCHECK, LLFOUND
INTEGER :: IPRINT 
 
LLDEFAULTVAL=.TRUE. 
LLREADNAM=.TRUE. 
LLCHECK=.TRUE. 
IPRINT=0 
IF(PRESENT(LDDEFAULTVAL)) LLDEFAULTVAL=LDDEFAULTVAL 
IF(PRESENT(LDREADNAM   )) LLREADNAM   =LDREADNAM 
IF(PRESENT(LDCHECK     )) LLCHECK     =LDCHECK 
IF(PRESENT(KPRINT      )) IPRINT      =KPRINT 
! 
!*      1. DEFAULT VALUES 
!       ----------------- 
! 
IF(LLDEFAULTVAL) THEN 
  !NOTES ON GENERAL DEFAULTS AND MODEL-SPECIFIC DEFAULTS :
  !- General default values *MUST* remain unchanged.
  !- To change the default value for a given application,                                                 
  !  an "IF(HPROGRAM=='...')" condition must be used.

  LNUCL=.TRUE.
  LSEDI=.TRUE.
  LHHONI = .FALSE.
  LMEYERS = .FALSE.
  NMOM_I = 2
  NMOM_S = 1
  NMOM_G = 1
  NMOM_H = 0
  NMOD_IFN = 1
  XIFN_CONC(:) = 100.
  LIFN_HOM = .TRUE.
  CIFN_SPECIES = 'PHILLIPS'
  CINT_MIXING = 'DM2'
  NMOD_IMM = 0
  NIND_SPECIE = 1
  LSNOW_T = .FALSE.
  CPRISTINE_ICE_LIMA = 'PLAT'
  CHEVRIMED_ICE_LIMA = 'GRAU'
  !XALPHAI=
  !XNUI=
  !XALPHAS=
  !XNUS=
  !XALPHAG=
  !XNUG=
  XFACTNUC_DEP = 1.0
  XFACTNUC_CON = 1.0
  NPHILLIPS=8
  LCIBU = .FALSE.
  XNDEBRIS_CIBU = 50.0
  LRDSF = .FALSE.
  LMURAKAMI=.TRUE.
  LACTI  = .TRUE.
  LSEDC  = .TRUE.
  LACTIT = .FALSE.
  LSPRO = .FALSE.
  LADJ   = .TRUE.
  LKHKO  = .FALSE.
  LKESSLERAC = .FALSE.
  NMOM_C = 2
  NMOM_R = 2
  NMOD_CCN = 1
  XCCN_CONC(:)=300.
  LCCN_HOM = .TRUE.
  CCCN_MODES = 'COPT'
  HINI_CCN   = 'AER'
  HTYPE_CCN(:) = 'M'
  XALPHAC = 3.0
  XNUC    = 1.0
  XALPHAR = 1.0
  XNUR    = 2.0
  XFSOLUB_CCN = 1.0
  XACTEMP_CCN = 280.
  XAERDIFF    = 0.0
  XAERHEIGHT  = 2000.
  LSCAV      = .FALSE.
  LAERO_MASS = .FALSE.
  LDEPOC = .TRUE.
  XVDEPOC = 0.02 ! 2 cm/s
  LACTTKE = .TRUE.
  LPTSPLIT     = .TRUE.
  LFEEDBACKT = .TRUE.
  NMAXITER  =  5
  XMRSTEP    = 0.005
  XTSTEP_TS  = 20.
ENDIF
!
!*      2. NAMELIST
!       -----------
!
IF(LLREADNAM) THEN
  CALL POSNAM_PHY(TFILENAM, 'NAM_PARAM_LIMA', LDNEEDNAM, LLFOUND)
  IF(LLFOUND) READ(UNIT=TFILENAM%NLU, NML=NAM_PARAM_LIMA)
ENDIF
!
!*      3. CHECKS
!       ---------
!
IF(LLCHECK) THEN
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CPRISTINE_ICE_LIMA', CPRISTINE_ICE_LIMA, &
                                                'PLAT', 'COLU', 'BURO')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CHEVRIMED_ICE_LIMA', CHEVRIMED_ICE_LIMA, &
                                                'GRAU', 'HAIL')

  IF ((LACTI .AND. HINI_CCN  == 'XXX')) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_LIMA', &
           &"YOU WANT TO USE A 2-MOMENT MICROPHYSICAL " // &
           &" SCHEME BUT YOU DIDNT FILL CORRECTLY NAM_PARAM_LIMA" // &
           &" YOU HAVE TO FILL HINI_CCN ")
  END IF

  IF(LACTI .AND. NMOD_CCN == 0) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_LIMA', &
           &"ACTIVATION OF AEROSOL PARTICLES IS NOT " // &
           &"POSSIBLE IF NMOD_CCN HAS VALUE ZERO. YOU HAVE TO SET AN UPPER " // &
           &"VALUE OF NMOD_CCN IN ORDER TO USE LIMA WARM ACTIVATION SCHEME.") 
  END IF

  IF(LNUCL .AND. NMOD_IFN == 0 .AND. (.NOT.LMEYERS) .AND. NMOM_I >= 2) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_LIMA', &
           &"NUCLEATION BY DEPOSITION AND CONTACT IS NOT " // &
           &"POSSIBLE IF NMOD_IFN HAS VALUE ZERO. YOU HAVE TO SET AN UPPER" //  &
           &"VALUE OF NMOD_IFN IN ORDER TO USE LIMA COLD NUCLEATION SCHEME.") 
  END IF

  IF(HPROGRAM=='AROME' .OR. HPROGRAM=='PHYEX') THEN
    IF(.NOT. LPTSPLIT) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_LIMA', &
           &"LPTSPLIT must be .TRUE. with this program: " // HPROGRAM)
    ENDIF
    IF(LSPRO) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_LIMA', &
           &"LSPRO must be .FALSE. with this program: " // HPROGRAM)
    ENDIF
  ENDIF
ENDIF
!
!*      3. PRINTS
!       ---------
!
IF(IPRINT>=1) THEN
  WRITE(UNIT=KLUOUT, NML=NAM_PARAM_LIMA)
ENDIF
!
END SUBROUTINE PARAM_LIMA_INIT
!
END MODULE MODD_PARAM_LIMA
