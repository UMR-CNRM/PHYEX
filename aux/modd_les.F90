!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODD_LES
!     ###############
!
!!****  *MODD_LES* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the
!     resolved fluxes and the spectra computed in LES mode
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_LES)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!
!!    AUTHOR
!!    ------
!!      J. Cuxart   *INM and Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    March 10, 1995
!!
!!       (J.Stein)  Sept. 25, 1995  add the model number in LES mode
!!       J. Cuxart  Oct.   4, 1996  New time series
!!       V. Masson  Jan.  20, 2000  New LES routines variables & //
!!       V. Masson  Nov.   6, 2002  LES budgets
!!       F. Couvreux Oct   1, 2006  LES PDF
!!       J.Pergaud   Oct    , 2007  MF LES
!!       P. Aumond   Oct     ,2009  User multimaskS + 4th order
!!       C.Lac       Oct     ,2014  Correction on user masks
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 30/03/2021: budgets: LES cartesian subdomain limits are defined in the physical domain
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
PUBLIC :: LES_ALLOCATE_DIM
INTERFACE LES_ALLOCATE_DIM
  MODULE PROCEDURE LES_ALLOCATE_1DIMX, LES_ALLOCATE_2DIMX, &
                   LES_ALLOCATE_3DIMX, LES_ALLOCATE_4DIMX, &
                   LES_ALLOCATE_3DIML, LES_ALLOCATE_4DIML, &
                   LES_ALLOCATE_3DIMI, LES_ALLOCATE_1DIMI,  &
                   LES_ALLOCATE_2DIMC
END INTERFACE LES_ALLOCATE_DIM

TYPE TLES_t
!-------------------------------------------------------------------------------
!
!* namelist variables
!
LOGICAL :: LLES_MEAN      ! flag to activate the mean computations
LOGICAL :: LLES_RESOLVED  ! flag to activate the resolved var. computations
LOGICAL :: LLES_SUBGRID   ! flag to activate the subgrid var.  computations
LOGICAL :: LLES_UPDRAFT   ! flag to activate the computations in updrafts
LOGICAL :: LLES_DOWNDRAFT ! flag to activate the computations in downdrafts
LOGICAL :: LLES_SPECTRA   ! flag to activate the spectra computations
LOGICAL :: LLES_PDF      ! flag to activate the pdf computations
!
INTEGER, DIMENSION(900) :: NLES_LEVELS         ! physical model levels for LES comp.
REAL,    DIMENSION(900) :: XLES_ALTITUDES      ! alt.  levels for LES comp.
INTEGER, DIMENSION(900) :: NSPECTRA_LEVELS     ! physical model levels for spectra comp.
REAL,    DIMENSION(900) :: XSPECTRA_ALTITUDES  ! alt.  levels for spectra comp.
!
INTEGER, DIMENSION(10) :: NLES_TEMP_SERIE_I   ! I, J and Z point
INTEGER, DIMENSION(10) :: NLES_TEMP_SERIE_J   ! localizations to
INTEGER, DIMENSION(10) :: NLES_TEMP_SERIE_Z   ! record temporal data

CHARACTER(LEN=4) :: CLES_NORM_TYPE ! type of turbulence normalization
CHARACTER(LEN=3) :: CBL_HEIGHT_DEF ! definition of the boundary layer height

REAL :: XLES_TEMP_SAMPLING    ! temporal sampling between each computation
REAL :: XLES_TEMP_MEAN_START  ! time (in s) from the beginning of the simulation
REAL :: XLES_TEMP_MEAN_END    ! for start and end of the temporal averaged comp.
REAL :: XLES_TEMP_MEAN_STEP   ! time step for each averaging

LOGICAL :: LLES_CART_MASK     ! flag to use a cartesian mask
INTEGER :: NLES_IINF          ! definition of the cartesians mask in physical domain
INTEGER :: NLES_ISUP          !     for NLES_CART_MODNBR model
INTEGER :: NLES_JINF          !               "
INTEGER :: NLES_JSUP          !               "
LOGICAL :: LLES_NEB_MASK      ! flag to use a 2D nebulosity mask
LOGICAL :: LLES_CORE_MASK     ! flag to use a 3D cloud core mask
LOGICAL :: LLES_MY_MASK       ! flag to use its own mask (must be coded by user)
INTEGER :: NLES_MASKS_USER    ! number of user masks for LES computations
LOGICAL :: LLES_CS_MASK       ! flag to use conditional sampling mask
INTEGER :: NPDF         ! number of pdf intervals
!
!-------------------------------------------------------------------------------
!
INTEGER, DIMENSION(JPMODELMAX) :: NLESn_IINF ! definition of the cartesians mask in physical domain
INTEGER, DIMENSION(JPMODELMAX) :: NLESn_ISUP !          for all models
INTEGER, DIMENSION(JPMODELMAX) :: NLESn_JINF !               "
INTEGER, DIMENSION(JPMODELMAX) :: NLESn_JSUP !               "
!
CHARACTER(LEN=4), DIMENSION(2,JPMODELMAX) :: CLES_LBCX
! X boundary conditions for 2 points correlations computations for all models
!
CHARACTER(LEN=4), DIMENSION(2,JPMODELMAX) :: CLES_LBCY
! Y boundary conditions for 2 points correlations computations for all models
!
!-------------------------------------------------------------------------------
!
LOGICAL :: LLES               ! flag to compute the LES diagnostics
!
LOGICAL :: LLES_CALL          ! flag to compute the LES diagnostics at current
!                             ! time step
!
!
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_CART_MASK
! 2D cartesian mask of the current model
!
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_NEB_MASK
! 2D nebulosity mask of the current model
!
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_CORE_MASK
! 2D surface precipitations mask of the current model
!
! 2D owner mask of the current model
LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: LLES_CURRENT_MY_MASKS
!
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_CS1_MASK
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_CS2_MASK
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LLES_CURRENT_CS3_MASK
! 2D conditional sampling mask of the current model
!
INTEGER :: NLES_CURRENT_TCOUNT
! current model LES time counter
!
INTEGER :: NLES_CURRENT_TIMES
! current model NLES_TIMES (number of LES samplings)
!
INTEGER :: NLES_CURRENT_IINF, NLES_CURRENT_ISUP, NLES_CURRENT_JINF, NLES_CURRENT_JSUP
! coordinates (in physical domain) for write_diachro, set to NLESn_IINF(current model), etc...
!
REAL :: XLES_CURRENT_DOMEGAX, XLES_CURRENT_DOMEGAY
! minimum wavelength in spectra analysis
!
CHARACTER(LEN=4), DIMENSION(2) :: CLES_CURRENT_LBCX
! current model X boundary conditions for 2 points correlations computations
!
CHARACTER(LEN=4), DIMENSION(2) :: CLES_CURRENT_LBCY
! current model Y boundary conditions for 2 points correlations computations
!
REAL, DIMENSION(:),   ALLOCATABLE :: XLES_CURRENT_Z
! altitudes for diachro
!
REAL :: XLES_CURRENT_ZS
! orography (used for normalization of altitudes)
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NKLIN_CURRENT_LES
! levels for vertical interpolation
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCOEFLIN_CURRENT_LES
! coefficients for vertical interpolation
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NKLIN_CURRENT_SPEC
! levels for vertical interpolation
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCOEFLIN_CURRENT_SPEC
! coefficients for vertical interpolation
!
REAL,DIMENSION(2) :: XTIME_LES
! time spent in subgrid LES computations in this time-step in TURB
!
!-------------------------------------------------------------------------------
!
!* normalization variables
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_M
! normalization coefficient for distances (Meters)
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_K
! normalization coefficient for temperatures (Kelvin)
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_S
! normalization coefficient for times (Seconds)
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_RHO
! normalization coefficient for densities
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_RV
! normalization coefficient for mixing ratio
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XLES_NORM_SV
! normalization coefficient for scalar variables
!
REAL, DIMENSION(:), ALLOCATABLE :: XLES_NORM_P
! normalization coefficient for pressure
!
!-------------------------------------------------------------------------------
!
!* monitoring variables
!
INTEGER :: NLES_MASKS         ! number of masks for LES computations
INTEGER :: NLES_K             ! number of vertical levels for local diagnostics
INTEGER :: NSPECTRA_K         ! number of vertical levels for spectra
!
CHARACTER(LEN=1) :: CLES_LEVEL_TYPE     ! type of vertical levels for local diag.
CHARACTER(LEN=1) :: CSPECTRA_LEVEL_TYPE ! type of vertical levels for spectra
!
!-------------------------------------------------------------------------------
!
!* subgrid variables for current model
!
!                                                                ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_W_SBG_WThl !  <w'w'Thl'>
!                                                                _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_W_SBG_WRt  !  <w'w'Rt'>
!                                                                _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_W_SBG_Thl2 !  <w'Thl'2>
!                                                                ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_W_SBG_Rt2  !  <w'Rt'2>
!                                                                _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_W_SBG_ThlRt!  <w'Thl'Rt'>
!                                                                _____
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_RES_W_SBG_WSv  !  <w'w'Sv'>
!                                                                ____
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_RES_W_SBG_Sv2 !  <w'Sv'2>
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XLES_SUBGRID_RCSIGS ! rc sigmas
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XLES_SUBGRID_RCSIGC ! rc sigmac
!                                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_U_SBG_UaU   !  <du'/dxa ua'u'>
!                                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_V_SBG_UaV   !  <dv'/dxa ua'v'>
!                                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_W_SBG_UaW   !  <dw'/dxa ua'w'>
!                                                                            _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_W_SBG_UaThl !  <dw'/dxa ua'Thl'>
!                                                                              _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Thl_SBG_UaW !  <dThl'/dxa ua'w'>
!                                                                              ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddz_Thl_SBG_W2   !  <dThl'/dz w'2>
!                                                                            ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_W_SBG_UaRt  !  <dw'/dxa ua'Rt'>
!                                                                             _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Rt_SBG_UaW  !  <dRt'/dxa ua'w'>
!                                                                             ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddz_Rt_SBG_W2    !  <dRt'/dz w'2>
!                                                                              ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Thl_SBG_UaRt!  <dThl'/dxa ua'Rt'>
!                                                                             _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Rt_SBG_UaThl!  <dRt'/dxa ua'Thl'>
!                                                                              _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Thl_SBG_UaThl! <dThl'/dxa ua'Thl'>
!                                                                             ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Rt_SBG_UaRt !  <dRt'/dxa ua'Rt'>
!                                                                              ______
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_W_SBG_UaSv  ! <dw'/dxa ua'Sv'>
!                                                                               _____
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Sv_SBG_UaW  ! <dSv'/dxa ua'w'>
!                                                                              ___
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_RES_ddz_Sv_SBG_W2    ! <dSv'/dz w'2>
!                                                                               ______
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_RES_ddxa_Sv_SBG_UaSv ! <dSv'/dxa ua'Sv'>
!
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_U2    ! <u'2>
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_V2    ! <v'2>
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_W2    ! <w'2>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Thl2  ! <Thl'2>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Rt2   ! <Rt'2>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Rc2   ! <Rc'2>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Ri2   ! <Ri'2>
!                                                            _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_ThlRt ! <Thl'Rt'>
!                                                             ____
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_Sv2   ! <Sv'2>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_UV    ! <u'v'>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WU    ! <w'u'>
!                                                            ____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WV    ! <w'v'>
!                                                            ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_UThl  ! <u'Thl'>
!                                                            ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_VThl  ! <v'Thl'>
!                                                            ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WThl  ! <w'Thl'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_URt   ! <u'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_VRt   ! <v'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WRt   ! <w'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_URc   ! <u'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_VRc   ! <v'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WRc   ! <w'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_SUBGRID_USv ! <u'Sv'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_SUBGRID_VSv ! <v'Sv'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WSv ! <w'Sv'>
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_UTke  ! <u'e>
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_VTke  ! <v'e>
!                                                            ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WTke  ! <w'e>
!                                                                   ___
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_ddz_WTke  !  <dw'e/dz>
!                                                             ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WThv   ! <w'Thv'>
!                                                             ________
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_ThlThv ! <Thl'Thv'>
!                                                             _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RtThv  ! <Rt'Thv'>
!                                                             _______
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_SvThv  ! <Sv'Thv'>
!                                                             ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_W2Thl  ! <w'2Thl>
!                                                              _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_W2Rt   ! <w'2Rt>
!                                                              _____
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_W2Sv   ! <w'2Sv>
!                                                             _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WThlRt ! <w'ThlRt>
!                                                             ______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WThl2  ! <w'Thl2>
!                                                             _____
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WRt2   ! <w'Rt2>
!                                                              _____
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_WSv2   ! <w'Sv2>
!                                                               _______
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_DISS_Tke ! <epsilon>
!                                                                 ____________
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_DISS_Thl2 ! <epsilon_Thl2>
!                                                                 ___________
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_DISS_Rt2  ! <epsilon_Rt2>
!                                                                ______________
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_DISS_ThlRt! <epsilon_ThlRt>
!                                                                ___________
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_DISS_Sv2  ! <epsilon_Sv2>
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WP     ! <w'p'>
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_ThlPz  ! <Thl'dp'/dz>
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RtPz   ! <Rt'dp'/dz>
!
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: X_LES_SUBGRID_SvPz   ! <Sv'dp'/dz>
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_PHI3   ! phi3
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_PSI3   ! psi3
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_LMix   ! mixing length
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_LDiss  ! dissipative length
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Km     ! eddy diffusivity for momentum
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_Kh     ! eddy diffusivity for heat
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_THLUP_MF  ! Thl of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RTUP_MF   ! Rt of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RVUP_MF   ! Rv of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RCUP_MF   ! Rc of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_RIUP_MF   ! Ri of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WUP_MF    ! Thl of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_MASSFLUX  ! Mass Flux
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_DETR      ! Detrainment
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_ENTR      ! Entrainment
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_FRACUP    ! Updraft Fraction
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_THVUP_MF  ! Thv of the Updraft
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WTHLMF ! Flux of thl
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WRTMF  ! Flux of rt
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WTHVMF ! Flux of thv
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WUMF   ! Flux of u
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_SUBGRID_WVMF   ! Flux of v
!
!* surface variables
!
REAL, DIMENSION(:), ALLOCATABLE :: X_LES_USTAR     ! local u* temporal series
REAL, DIMENSION(:), ALLOCATABLE :: X_LES_UW0       ! uw temporal series
REAL, DIMENSION(:), ALLOCATABLE :: X_LES_VW0       ! vw temporal series
REAL, DIMENSION(:), ALLOCATABLE :: X_LES_Q0        ! Qo temporal series
REAL, DIMENSION(:), ALLOCATABLE :: X_LES_E0        ! Eo temporal series
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_SV0     ! scalar surface fluxes
!
!* pdf variables
REAL                              :: XRV_PDF_MIN       ! min of rv pdf
REAL                              :: XRV_PDF_MAX       ! max of rv pdf
REAL                              :: XTH_PDF_MIN       ! min of theta pdf
REAL                              :: XTH_PDF_MAX       ! max of theta pdf
REAL                              :: XW_PDF_MIN       ! min of w pdf
REAL                              :: XW_PDF_MAX       ! max of w pdf
REAL                              :: XTHV_PDF_MIN       ! min of thetav pdf
REAL                              :: XTHV_PDF_MAX       ! max of thetav pdf
REAL                              :: XRC_PDF_MIN       ! min of rc pdf
REAL                              :: XRC_PDF_MAX       ! max of rc pdf
REAL                              :: XRR_PDF_MIN       ! min of rr pdf
REAL                              :: XRR_PDF_MAX       ! max of rr pdf
REAL                              :: XRI_PDF_MIN       ! min of ri pdf
REAL                              :: XRI_PDF_MAX       ! max of ri pdf
REAL                              :: XRS_PDF_MIN       ! min of rs pdf
REAL                              :: XRS_PDF_MAX       ! max of rs pdf
REAL                              :: XRG_PDF_MIN       ! min of rg pdf
REAL                              :: XRG_PDF_MAX       ! max of rg pdf
REAL                              :: XRT_PDF_MIN       ! min of rt pdf
REAL                              :: XRT_PDF_MAX       ! max of rt pdf
REAL                              :: XTHL_PDF_MIN       ! min of thetal pdf
REAL                              :: XTHL_PDF_MAX       ! max of thetal pdf
!-------------------------------------------------------------------------------
!* pdf distribution
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RV   ! rv  pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_TH   ! theta pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_W    ! w  pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_THV  ! thetav pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RC   ! rc  pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RR   ! rr pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RI   ! ri  pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RS   ! rs pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RG   ! rg pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_RT   ! rt  pdf
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XLES_PDF_THL  ! thetal pdf
!
END TYPE TLES_t
!
TYPE(TLES_t), SAVE, TARGET :: TLES
!
!-------------------------------------------------------------------------------
!
!* namelist variables
!
LOGICAL, POINTER :: LLES_MEAN     => NULL() ! flag to activate the mean computations
LOGICAL, POINTER :: LLES_RESOLVED => NULL() ! flag to activate the resolved var. computations
LOGICAL, POINTER :: LLES_SUBGRID  => NULL() ! flag to activate the subgrid var.  computations
LOGICAL, POINTER :: LLES_UPDRAFT  => NULL() ! flag to activate the computations in updrafts
LOGICAL, POINTER :: LLES_DOWNDRAFT=> NULL() ! flag to activate the computations in downdrafts
LOGICAL, POINTER :: LLES_SPECTRA  => NULL() ! flag to activate the spectra computations
LOGICAL, POINTER :: LLES_PDF     => NULL() ! flag to activate the pdf computations
!
INTEGER, DIMENSION(:), POINTER :: NLES_LEVELS        => NULL() ! physical model levels for LES comp.
REAL,    DIMENSION(:), POINTER :: XLES_ALTITUDES     => NULL() ! alt.  levels for LES comp.
INTEGER, DIMENSION(:), POINTER :: NSPECTRA_LEVELS    => NULL() ! physical model levels for spectra comp.
REAL,    DIMENSION(:), POINTER :: XSPECTRA_ALTITUDES => NULL() ! alt.  levels for spectra comp.
!
INTEGER, DIMENSION(:), POINTER :: NLES_TEMP_SERIE_I  => NULL() ! I, J and Z point
INTEGER, DIMENSION(:), POINTER :: NLES_TEMP_SERIE_J  => NULL() ! localizations to
INTEGER, DIMENSION(:), POINTER :: NLES_TEMP_SERIE_Z  => NULL() ! record temporal data

CHARACTER(LEN=4), POINTER :: CLES_NORM_TYPE=> NULL() ! type of turbulence normalization
CHARACTER(LEN=3), POINTER :: CBL_HEIGHT_DEF=> NULL() ! definition of the boundary layer height

REAL, POINTER :: XLES_TEMP_SAMPLING   => NULL() ! temporal sampling between each computation
REAL, POINTER :: XLES_TEMP_MEAN_START => NULL() ! time (in s) from the beginning of the simulation
REAL, POINTER :: XLES_TEMP_MEAN_END   => NULL() ! for start and end of the temporal averaged comp.
REAL, POINTER :: XLES_TEMP_MEAN_STEP  => NULL() ! time step for each averaging

LOGICAL, POINTER :: LLES_CART_MASK    => NULL() ! flag to use a cartesian mask
INTEGER, POINTER :: NLES_IINF         => NULL() ! definition of the cartesians mask in physical domain
INTEGER, POINTER :: NLES_ISUP         => NULL() !     for NLES_CART_MODNBR model
INTEGER, POINTER :: NLES_JINF         => NULL() !               "
INTEGER, POINTER :: NLES_JSUP         => NULL() !               "
LOGICAL, POINTER :: LLES_NEB_MASK     => NULL() ! flag to use a 2D nebulosity mask
LOGICAL, POINTER :: LLES_CORE_MASK    => NULL() ! flag to use a 3D cloud core mask
LOGICAL, POINTER :: LLES_MY_MASK      => NULL() ! flag to use its own mask (must be coded by user)
INTEGER, POINTER :: NLES_MASKS_USER   => NULL() ! number of user masks for LES computations
LOGICAL, POINTER :: LLES_CS_MASK      => NULL() ! flag to use conditional sampling mask
INTEGER, POINTER :: NPDF        => NULL() ! number of pdf intervals
!
!-------------------------------------------------------------------------------
!
INTEGER, DIMENSION(:), POINTER :: NLESn_IINF=> NULL() ! definition of the cartesians mask in physical domain
INTEGER, DIMENSION(:), POINTER :: NLESn_ISUP=> NULL() !          for all models
INTEGER, DIMENSION(:), POINTER :: NLESn_JINF=> NULL() !               "
INTEGER, DIMENSION(:), POINTER :: NLESn_JSUP=> NULL() !               "
!
CHARACTER(LEN=4), DIMENSION(:,:), POINTER :: CLES_LBCX=> NULL()
! X boundary conditions for 2 points correlations computations for all models
!
CHARACTER(LEN=4), DIMENSION(:,:), POINTER :: CLES_LBCY=> NULL()
! Y boundary conditions for 2 points correlations computations for all models
!
!-------------------------------------------------------------------------------
!
LOGICAL, POINTER :: LLES              => NULL() ! flag to compute the LES diagnostics
!
LOGICAL, POINTER :: LLES_CALL         => NULL() ! flag to compute the LES diagnostics at current
!                            => NULL() ! time step
!
!
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_CART_MASK=> NULL()
! 2D cartesian mask of the current model
!
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_NEB_MASK=> NULL()
! 2D nebulosity mask of the current model
!
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_CORE_MASK=> NULL()
! 2D surface precipitations mask of the current model
!
! 2D owner mask of the current model
LOGICAL, DIMENSION(:,:,:,:), POINTER :: LLES_CURRENT_MY_MASKS=> NULL()
!
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_CS1_MASK=> NULL()
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_CS2_MASK=> NULL()
LOGICAL, DIMENSION(:,:,:), POINTER :: LLES_CURRENT_CS3_MASK=> NULL()
! 2D conditional sampling mask of the current model
!
INTEGER, POINTER :: NLES_CURRENT_TCOUNT=> NULL()
! current model LES time counter
!
INTEGER, POINTER :: NLES_CURRENT_TIMES=> NULL()
! current model NLES_TIMES (number of LES samplings)
!
INTEGER, POINTER :: NLES_CURRENT_IINF=> NULL(), NLES_CURRENT_ISUP=> NULL(), &
                    NLES_CURRENT_JINF=> NULL(), NLES_CURRENT_JSUP=> NULL()
! coordinates (in physical domain) for write_diachro, set to NLESn_IINF(current model), etc...
!
REAL, POINTER :: XLES_CURRENT_DOMEGAX=> NULL(), XLES_CURRENT_DOMEGAY=> NULL()
! minimum wavelength in spectra analysis
!
CHARACTER(LEN=4), DIMENSION(:), POINTER :: CLES_CURRENT_LBCX=> NULL()
! current model X boundary conditions for 2 points correlations computations
!
CHARACTER(LEN=4), DIMENSION(:), POINTER :: CLES_CURRENT_LBCY=> NULL()
! current model Y boundary conditions for 2 points correlations computations
!
REAL, DIMENSION(:),   POINTER :: XLES_CURRENT_Z=> NULL()
! altitudes for diachro
!
REAL, POINTER :: XLES_CURRENT_ZS=> NULL()
! orography (used for normalization of altitudes)
!
INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_CURRENT_LES=> NULL()
! levels for vertical interpolation
!
REAL, DIMENSION(:,:,:), POINTER :: XCOEFLIN_CURRENT_LES=> NULL()
! coefficients for vertical interpolation
!
INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_CURRENT_SPEC=> NULL()
! levels for vertical interpolation
!
REAL, DIMENSION(:,:,:), POINTER :: XCOEFLIN_CURRENT_SPEC=> NULL()
! coefficients for vertical interpolation
!
REAL,DIMENSION(:), POINTER :: XTIME_LES=> NULL()
! time spent in subgrid LES computations in this time-step in TURB
!
!-------------------------------------------------------------------------------
!
!* normalization variables
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_M=> NULL()
! normalization coefficient for distances (Meters)
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_K=> NULL()
! normalization coefficient for temperatures (Kelvin)
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_S=> NULL()
! normalization coefficient for times (Seconds)
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_RHO=> NULL()
! normalization coefficient for densities
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_RV=> NULL()
! normalization coefficient for mixing ratio
!
REAL, DIMENSION(:,:), POINTER :: XLES_NORM_SV=> NULL()
! normalization coefficient for scalar variables
!
REAL, DIMENSION(:), POINTER :: XLES_NORM_P=> NULL()
! normalization coefficient for pressure
!
!-------------------------------------------------------------------------------
!
!* monitoring variables
!
INTEGER, POINTER :: NLES_MASKS        => NULL() ! number of masks for LES computations
INTEGER, POINTER :: NLES_K            => NULL() ! number of vertical levels for local diagnostics
INTEGER, POINTER :: NSPECTRA_K        => NULL() ! number of vertical levels for spectra
!
CHARACTER(LEN=1), POINTER :: CLES_LEVEL_TYPE    => NULL() ! type of vertical levels for local diag.
CHARACTER(LEN=1), POINTER :: CSPECTRA_LEVEL_TYPE=> NULL() ! type of vertical levels for spectra
!
!-------------------------------------------------------------------------------
!
!* subgrid variables for current model
!
!                                                                ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_W_SBG_WThl=> NULL() !  <w'w'Thl'>
!                                                                _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_W_SBG_WRt => NULL() !  <w'w'Rt'>
!                                                                _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_W_SBG_Thl2=> NULL() !  <w'Thl'2>
!                                                                ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_W_SBG_Rt2 => NULL() !  <w'Rt'2>
!                                                                _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_W_SBG_ThlRt=> NULL()!  <w'Thl'Rt'>
!                                                                _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_W_SBG_WSv => NULL() !  <w'w'Sv'>
!                                                                ____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_W_SBG_Sv2=> NULL() !  <w'Sv'2>
!
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RCSIGS=> NULL() ! rc sigmas
!
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RCSIGC=> NULL() ! rc sigmac
!                                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_U_SBG_UaU  => NULL() !  <du'/dxa ua'u'>
!                                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_V_SBG_UaV  => NULL() !  <dv'/dxa ua'v'>
!                                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_W_SBG_UaW  => NULL() !  <dw'/dxa ua'w'>
!                                                                            _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_W_SBG_UaThl=> NULL() !  <dw'/dxa ua'Thl'>
!                                                                              _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Thl_SBG_UaW=> NULL() !  <dThl'/dxa ua'w'>
!                                                                              ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddz_Thl_SBG_W2  => NULL() !  <dThl'/dz w'2>
!                                                                            ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_W_SBG_UaRt => NULL() !  <dw'/dxa ua'Rt'>
!                                                                             _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Rt_SBG_UaW => NULL() !  <dRt'/dxa ua'w'>
!                                                                             ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddz_Rt_SBG_W2   => NULL() !  <dRt'/dz w'2>
!                                                                              ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Thl_SBG_UaRt=> NULL()!  <dThl'/dxa ua'Rt'>
!                                                                             _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Rt_SBG_UaThl=> NULL()!  <dRt'/dxa ua'Thl'>
!                                                                              _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Thl_SBG_UaThl=> NULL()! <dThl'/dxa ua'Thl'>
!                                                                             ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_RES_ddxa_Rt_SBG_UaRt=> NULL() !  <dRt'/dxa ua'Rt'>
!                                                                              ______
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_ddxa_W_SBG_UaSv => NULL() ! <dw'/dxa ua'Sv'>
!                                                                               _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_ddxa_Sv_SBG_UaW => NULL() ! <dSv'/dxa ua'w'>
!                                                                              ___
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_ddz_Sv_SBG_W2   => NULL() ! <dSv'/dz w'2>
!                                                                               ______
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_RES_ddxa_Sv_SBG_UaSv=> NULL() ! <dSv'/dxa ua'Sv'>
!
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_U2   => NULL() ! <u'2>
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_V2   => NULL() ! <v'2>
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_W2   => NULL() ! <w'2>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Thl2 => NULL() ! <Thl'2>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Rt2  => NULL() ! <Rt'2>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Rc2  => NULL() ! <Rc'2>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Ri2  => NULL() ! <Ri'2>
!                                                            _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_ThlRt=> NULL() ! <Thl'Rt'>
!                                                             ____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_Sv2  => NULL() ! <Sv'2>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_UV   => NULL() ! <u'v'>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WU   => NULL() ! <w'u'>
!                                                            ____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WV   => NULL() ! <w'v'>
!                                                            ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_UThl => NULL() ! <u'Thl'>
!                                                            ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_VThl => NULL() ! <v'Thl'>
!                                                            ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WThl => NULL() ! <w'Thl'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_URt  => NULL() ! <u'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_VRt  => NULL() ! <v'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WRt  => NULL() ! <w'Rt'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_URc  => NULL() ! <u'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_VRc  => NULL() ! <v'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WRc  => NULL() ! <w'Rc'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_USv=> NULL() ! <u'Sv'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_VSv=> NULL() ! <v'Sv'>
!                                                            _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_WSv=> NULL() ! <w'Sv'>
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_UTke => NULL() ! <u'e>
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_VTke => NULL() ! <v'e>
!                                                            ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WTke => NULL() ! <w'e>
!                                                                   ___
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_ddz_WTke => NULL() !  <dw'e/dz>
!                                                             ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WThv  => NULL() ! <w'Thv'>
!                                                             ________
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_ThlThv=> NULL() ! <Thl'Thv'>
!                                                             _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RtThv => NULL() ! <Rt'Thv'>
!                                                             _______
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_SvThv => NULL() ! <Sv'Thv'>
!                                                             ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_W2Thl => NULL() ! <w'2Thl>
!                                                              _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_W2Rt  => NULL() ! <w'2Rt>
!                                                              _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_W2Sv  => NULL() ! <w'2Sv>
!                                                             _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WThlRt=> NULL() ! <w'ThlRt>
!                                                             ______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WThl2 => NULL() ! <w'Thl2>
!                                                             _____
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WRt2  => NULL() ! <w'Rt2>
!                                                              _____
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_WSv2  => NULL() ! <w'Sv2>
!                                                               _______
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_DISS_Tke=> NULL() ! <epsilon>
!                                                                 ____________
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_DISS_Thl2=> NULL() ! <epsilon_Thl2>
!                                                                 ___________
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_DISS_Rt2 => NULL() ! <epsilon_Rt2>
!                                                                ______________
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_DISS_ThlRt=> NULL()! <epsilon_ThlRt>
!                                                                ___________
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_DISS_Sv2 => NULL() ! <epsilon_Sv2>
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WP    => NULL() ! <w'p'>
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_ThlPz => NULL() ! <Thl'dp'/dz>
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RtPz  => NULL() ! <Rt'dp'/dz>
!
REAL, DIMENSION(:,:,:,:), POINTER :: X_LES_SUBGRID_SvPz  => NULL() ! <Sv'dp'/dz>
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_PHI3  => NULL() ! phi3
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_PSI3  => NULL() ! psi3
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_LMix  => NULL() ! mixing length
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_LDiss => NULL() ! dissipative length
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Km    => NULL() ! eddy diffusivity for momentum
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_Kh    => NULL() ! eddy diffusivity for heat
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_THLUP_MF => NULL() ! Thl of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RTUP_MF  => NULL() ! Rt of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RVUP_MF  => NULL() ! Rv of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RCUP_MF  => NULL() ! Rc of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_RIUP_MF  => NULL() ! Ri of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WUP_MF   => NULL() ! Thl of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_MASSFLUX => NULL() ! Mass Flux
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_DETR     => NULL() ! Detrainment
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_ENTR     => NULL() ! Entrainment
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_FRACUP   => NULL() ! Updraft Fraction
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_THVUP_MF => NULL() ! Thv of the Updraft
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WTHLMF=> NULL() ! Flux of thl
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WRTMF => NULL() ! Flux of rt
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WTHVMF=> NULL() ! Flux of thv
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WUMF  => NULL() ! Flux of u
!
REAL, DIMENSION(:,:,:), POINTER :: X_LES_SUBGRID_WVMF  => NULL() ! Flux of v
!
!* surface variables
!
REAL, DIMENSION(:), POINTER :: X_LES_USTAR    => NULL() ! local u* temporal series
REAL, DIMENSION(:), POINTER :: X_LES_UW0      => NULL() ! uw temporal series
REAL, DIMENSION(:), POINTER :: X_LES_VW0      => NULL() ! vw temporal series
REAL, DIMENSION(:), POINTER :: X_LES_Q0       => NULL() ! Qo temporal series
REAL, DIMENSION(:), POINTER :: X_LES_E0       => NULL() ! Eo temporal series
REAL, DIMENSION(:,:), POINTER :: X_LES_SV0    => NULL() ! scalar surface fluxes
!
!* pdf variables
REAL                             , POINTER :: XRV_PDF_MIN      => NULL() ! min of rv pdf
REAL                             , POINTER :: XRV_PDF_MAX      => NULL() ! max of rv pdf
REAL                             , POINTER :: XTH_PDF_MIN      => NULL() ! min of theta pdf
REAL                             , POINTER :: XTH_PDF_MAX      => NULL() ! max of theta pdf
REAL                             , POINTER :: XW_PDF_MIN       => NULL() ! min of w pdf
REAL                             , POINTER :: XW_PDF_MAX       => NULL() ! max of w pdf
REAL                             , POINTER :: XTHV_PDF_MIN     => NULL() ! min of thetav pdf
REAL                             , POINTER :: XTHV_PDF_MAX     => NULL() ! max of thetav pdf
REAL                             , POINTER :: XRC_PDF_MIN      => NULL() ! min of rc pdf
REAL                             , POINTER :: XRC_PDF_MAX      => NULL() ! max of rc pdf
REAL                             , POINTER :: XRR_PDF_MIN      => NULL() ! min of rr pdf
REAL                             , POINTER :: XRR_PDF_MAX      => NULL() ! max of rr pdf
REAL                             , POINTER :: XRI_PDF_MIN      => NULL() ! min of ri pdf
REAL                             , POINTER :: XRI_PDF_MAX      => NULL() ! max of ri pdf
REAL                             , POINTER :: XRS_PDF_MIN      => NULL() ! min of rs pdf
REAL                             , POINTER :: XRS_PDF_MAX      => NULL() ! max of rs pdf
REAL                             , POINTER :: XRG_PDF_MIN      => NULL() ! min of rg pdf
REAL                             , POINTER :: XRG_PDF_MAX      => NULL() ! max of rg pdf
REAL                             , POINTER :: XRT_PDF_MIN      => NULL() ! min of rt pdf
REAL                             , POINTER :: XRT_PDF_MAX      => NULL() ! max of rt pdf
REAL                             , POINTER :: XTHL_PDF_MIN     => NULL() ! min of thetal pdf
REAL                             , POINTER :: XTHL_PDF_MAX     => NULL() ! max of thetal pdf
!-------------------------------------------------------------------------------
!* pdf distribution
!
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RV  => NULL() ! rv  pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_TH  => NULL() ! theta pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_W   => NULL() ! w  pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_THV => NULL() ! thetav pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RC  => NULL() ! rc  pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RR  => NULL() ! rr pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RI  => NULL() ! ri  pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RS  => NULL() ! rs pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RG  => NULL() ! rg pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_RT  => NULL() ! rt  pdf
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_PDF_THL => NULL() ! thetal pdf
!-------------------------------------------------------------------------------
!!
CONTAINS
SUBROUTINE LES_ASSOCIATE()
  ! Associate all LES non-allocatable variables to the TYPE LES
  IMPLICIT NONE
  NLES_LEVELS => TLES%NLES_LEVELS
  XLES_ALTITUDES => TLES%XLES_ALTITUDES
  NSPECTRA_LEVELS => TLES%NSPECTRA_LEVELS
  XSPECTRA_ALTITUDES => TLES%XSPECTRA_ALTITUDES
  XTIME_LES => TLES%XTIME_LES
  CLES_LBCX => TLES%CLES_LBCX
  CLES_LBCY => TLES%CLES_LBCY
  CLES_CURRENT_LBCY => TLES%CLES_CURRENT_LBCY
  CLES_CURRENT_LBCX => TLES%CLES_CURRENT_LBCX
  NLESn_IINF => TLES%NLESn_IINF
  NLESn_ISUP => TLES%NLESn_ISUP
  NLESn_JINF => TLES%NLESn_JINF
  NLESn_JSUP => TLES%NLESn_JSUP
  NLES_TEMP_SERIE_I => TLES%NLES_TEMP_SERIE_I
  NLES_TEMP_SERIE_J => TLES%NLES_TEMP_SERIE_J
  NLES_TEMP_SERIE_Z => TLES%NLES_TEMP_SERIE_Z
  LLES_MEAN => TLES%LLES_MEAN
  LLES_RESOLVED => TLES%LLES_RESOLVED
  LLES_SUBGRID => TLES%LLES_SUBGRID
  LLES_UPDRAFT => TLES%LLES_UPDRAFT
  LLES_DOWNDRAFT => TLES%LLES_DOWNDRAFT
  LLES_SPECTRA => TLES%LLES_SPECTRA
  LLES_PDF => TLES%LLES_PDF
  CLES_NORM_TYPE => TLES%CLES_NORM_TYPE
  CBL_HEIGHT_DEF => TLES%CBL_HEIGHT_DEF
  XLES_TEMP_SAMPLING => TLES%XLES_TEMP_SAMPLING
  XLES_TEMP_MEAN_START => TLES%XLES_TEMP_MEAN_START
  XLES_TEMP_MEAN_END => TLES%XLES_TEMP_MEAN_END
  XLES_TEMP_MEAN_STEP => TLES%XLES_TEMP_MEAN_STEP
  LLES_CART_MASK => TLES%LLES_CART_MASK
  NLES_IINF => TLES%NLES_IINF
  NLES_ISUP => TLES%NLES_ISUP
  NLES_JINF => TLES%NLES_JINF
  NLES_JSUP => TLES%NLES_JSUP
  LLES_NEB_MASK => TLES%LLES_NEB_MASK
  LLES_CORE_MASK => TLES%LLES_CORE_MASK
  LLES_MY_MASK => TLES%LLES_MY_MASK
  NLES_MASKS_USER => TLES%NLES_MASKS_USER
  LLES_CS_MASK => TLES%LLES_CS_MASK
  NPDF => TLES%NPDF
  LLES => TLES%LLES
  LLES_CALL => TLES%LLES_CALL
  NLES_CURRENT_TCOUNT => TLES%NLES_CURRENT_TCOUNT
  NLES_CURRENT_TIMES => TLES%NLES_CURRENT_TIMES
  NLES_CURRENT_IINF => TLES%NLES_CURRENT_IINF
  NLES_CURRENT_ISUP => TLES%NLES_CURRENT_ISUP
  NLES_CURRENT_JINF => TLES%NLES_CURRENT_JINF
  NLES_CURRENT_JSUP => TLES%NLES_CURRENT_JSUP
  XLES_CURRENT_DOMEGAX => TLES%XLES_CURRENT_DOMEGAX
  XLES_CURRENT_DOMEGAY => TLES%XLES_CURRENT_DOMEGAY
  XLES_CURRENT_ZS => TLES%XLES_CURRENT_ZS
  NLES_MASKS => TLES%NLES_MASKS
  NLES_K => TLES%NLES_K
  NSPECTRA_K => TLES%NSPECTRA_K
  CLES_LEVEL_TYPE => TLES%CLES_LEVEL_TYPE
  CSPECTRA_LEVEL_TYPE => TLES%CSPECTRA_LEVEL_TYPE
  XRV_PDF_MIN => TLES%XRV_PDF_MIN
  XRV_PDF_MAX => TLES%XRV_PDF_MAX
  XTH_PDF_MIN => TLES%XTH_PDF_MIN
  XTH_PDF_MAX => TLES%XTH_PDF_MAX
  XW_PDF_MIN => TLES%XW_PDF_MIN
  XW_PDF_MAX => TLES%XW_PDF_MAX
  XTHV_PDF_MIN => TLES%XTHV_PDF_MIN
  XTHV_PDF_MAX => TLES%XTHV_PDF_MAX
  XRC_PDF_MIN => TLES%XRC_PDF_MIN
  XRC_PDF_MAX => TLES%XRC_PDF_MAX
  XRR_PDF_MIN => TLES%XRR_PDF_MIN
  XRR_PDF_MAX => TLES%XRR_PDF_MAX
  XRI_PDF_MIN => TLES%XRI_PDF_MIN
  XRI_PDF_MAX => TLES%XRI_PDF_MAX
  XRS_PDF_MIN => TLES%XRS_PDF_MIN
  XRS_PDF_MAX => TLES%XRS_PDF_MAX
  XRG_PDF_MIN => TLES%XRG_PDF_MIN
  XRG_PDF_MAX => TLES%XRG_PDF_MAX
  XRT_PDF_MIN => TLES%XRT_PDF_MIN
  XRT_PDF_MAX => TLES%XRT_PDF_MAX
  XTHL_PDF_MIN => TLES%XTHL_PDF_MIN
  XTHL_PDF_MAX => TLES%XTHL_PDF_MAX
END SUBROUTINE LES_ASSOCIATE
!
SUBROUTINE LES_ALLOCATE(HNAME,NDIMS)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, DIMENSION(:) :: NDIMS
  !
  SELECT CASE(HNAME)
  !
  CASE('LLES_CURRENT_CART_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_CART_MASK,NDIMS)
    LLES_CURRENT_CART_MASK=>TLES%LLES_CURRENT_CART_MASK
  CASE('LLES_CURRENT_NEB_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_NEB_MASK,NDIMS)
    LLES_CURRENT_NEB_MASK=>TLES%LLES_CURRENT_NEB_MASK
  CASE('LLES_CURRENT_CORE_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_CORE_MASK,NDIMS)
    LLES_CURRENT_CORE_MASK=>TLES%LLES_CURRENT_CORE_MASK
  CASE('LLES_CURRENT_MY_MASKS')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_MY_MASKS,NDIMS)
    LLES_CURRENT_MY_MASKS=>TLES%LLES_CURRENT_MY_MASKS
  CASE('LLES_CURRENT_CS1_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_CS1_MASK,NDIMS)
    LLES_CURRENT_CS1_MASK=>TLES%LLES_CURRENT_CS1_MASK
  CASE('LLES_CURRENT_CS2_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_CS2_MASK,NDIMS)
    LLES_CURRENT_CS2_MASK=>TLES%LLES_CURRENT_CS2_MASK
  CASE('LLES_CURRENT_CS3_MASK')
    CALL LES_ALLOCATE_DIM(TLES%LLES_CURRENT_CS3_MASK,NDIMS)
    LLES_CURRENT_CS3_MASK=>TLES%LLES_CURRENT_CS3_MASK
  CASE('XLES_CURRENT_Z')
    CALL LES_ALLOCATE_DIM(TLES%XLES_CURRENT_Z,NDIMS)
    XLES_CURRENT_Z=>TLES%XLES_CURRENT_Z
  CASE('NKLIN_CURRENT_LES')
    CALL LES_ALLOCATE_DIM(TLES%NKLIN_CURRENT_LES,NDIMS)
    NKLIN_CURRENT_LES=>TLES%NKLIN_CURRENT_LES
  CASE('XCOEFLIN_CURRENT_LES')
    CALL LES_ALLOCATE_DIM(TLES%XCOEFLIN_CURRENT_LES,NDIMS)
    XCOEFLIN_CURRENT_LES=>TLES%XCOEFLIN_CURRENT_LES
  CASE('NKLIN_CURRENT_SPEC')
    CALL LES_ALLOCATE_DIM(TLES%NKLIN_CURRENT_SPEC,NDIMS)
    NKLIN_CURRENT_SPEC=>TLES%NKLIN_CURRENT_SPEC
  CASE('XCOEFLIN_CURRENT_SPEC')
    CALL LES_ALLOCATE_DIM(TLES%XCOEFLIN_CURRENT_SPEC,NDIMS)
    XCOEFLIN_CURRENT_SPEC=>TLES%XCOEFLIN_CURRENT_SPEC
  CASE('XLES_NORM_M')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_M,NDIMS)
    XLES_NORM_M=>TLES%XLES_NORM_M
  CASE('XLES_NORM_K')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_K,NDIMS)
    XLES_NORM_K=>TLES%XLES_NORM_K
  CASE('XLES_NORM_S')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_S,NDIMS)
    XLES_NORM_S=>TLES%XLES_NORM_S
  CASE('XLES_NORM_RHO')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_RHO,NDIMS)
    XLES_NORM_RHO=>TLES%XLES_NORM_RHO
  CASE('XLES_NORM_RV')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_RV,NDIMS)
    XLES_NORM_RV=>TLES%XLES_NORM_RV
  CASE('XLES_NORM_SV')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_SV,NDIMS)
    XLES_NORM_SV=>TLES%XLES_NORM_SV
  CASE('XLES_NORM_P')
    CALL LES_ALLOCATE_DIM(TLES%XLES_NORM_P,NDIMS)
    XLES_NORM_P=>TLES%XLES_NORM_P
  CASE('X_LES_RES_W_SBG_WThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_WThl,NDIMS)
    X_LES_RES_W_SBG_WThl=>TLES%X_LES_RES_W_SBG_WThl
  CASE('X_LES_RES_W_SBG_WRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_WRt,NDIMS)
    X_LES_RES_W_SBG_WRt=>TLES%X_LES_RES_W_SBG_WRt
  CASE('X_LES_RES_W_SBG_Thl2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_Thl2,NDIMS)
    X_LES_RES_W_SBG_Thl2=>TLES%X_LES_RES_W_SBG_Thl2
  CASE('X_LES_RES_W_SBG_Rt2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_Rt2,NDIMS)
    X_LES_RES_W_SBG_Rt2=>TLES%X_LES_RES_W_SBG_Rt2
  CASE('X_LES_RES_W_SBG_ThlRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_ThlRt,NDIMS)
    X_LES_RES_W_SBG_ThlRt=>TLES%X_LES_RES_W_SBG_ThlRt
  CASE('X_LES_RES_W_SBG_WSv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_WSv,NDIMS)
    X_LES_RES_W_SBG_WSv=>TLES%X_LES_RES_W_SBG_WSv
  CASE('X_LES_RES_W_SBG_Sv2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_W_SBG_Sv2,NDIMS)
    X_LES_RES_W_SBG_Sv2=>TLES%X_LES_RES_W_SBG_Sv2
  CASE('XLES_SUBGRID_RCSIGS')
    CALL LES_ALLOCATE_DIM(TLES%XLES_SUBGRID_RCSIGS,NDIMS)
    XLES_SUBGRID_RCSIGS=>TLES%XLES_SUBGRID_RCSIGS
  CASE('XLES_SUBGRID_RCSIGC')
    CALL LES_ALLOCATE_DIM(TLES%XLES_SUBGRID_RCSIGC,NDIMS)
    XLES_SUBGRID_RCSIGC=>TLES%XLES_SUBGRID_RCSIGC
  CASE('X_LES_RES_ddxa_U_SBG_UaU')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_U_SBG_UaU,NDIMS)
    X_LES_RES_ddxa_U_SBG_UaU=>TLES%X_LES_RES_ddxa_U_SBG_UaU
  CASE('X_LES_RES_ddxa_V_SBG_UaV')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_V_SBG_UaV,NDIMS)
    X_LES_RES_ddxa_V_SBG_UaV=>TLES%X_LES_RES_ddxa_V_SBG_UaV
  CASE('X_LES_RES_ddxa_W_SBG_UaW')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_W_SBG_UaW,NDIMS)
    X_LES_RES_ddxa_W_SBG_UaW=>TLES%X_LES_RES_ddxa_W_SBG_UaW
  CASE('X_LES_RES_ddxa_W_SBG_UaThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_W_SBG_UaThl,NDIMS)
    X_LES_RES_ddxa_W_SBG_UaThl=>TLES%X_LES_RES_ddxa_W_SBG_UaThl
  CASE('X_LES_RES_ddxa_Thl_SBG_UaW')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Thl_SBG_UaW,NDIMS)
    X_LES_RES_ddxa_Thl_SBG_UaW=>TLES%X_LES_RES_ddxa_Thl_SBG_UaW
  CASE('X_LES_RES_ddz_Thl_SBG_W2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddz_Thl_SBG_W2,NDIMS)
    X_LES_RES_ddz_Thl_SBG_W2=>TLES%X_LES_RES_ddz_Thl_SBG_W2
  CASE('X_LES_RES_ddxa_W_SBG_UaRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_W_SBG_UaRt,NDIMS)
    X_LES_RES_ddxa_W_SBG_UaRt=>TLES%X_LES_RES_ddxa_W_SBG_UaRt
  CASE('X_LES_RES_ddxa_Rt_SBG_UaW')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Rt_SBG_UaW,NDIMS)
    X_LES_RES_ddxa_Rt_SBG_UaW=>TLES%X_LES_RES_ddxa_Rt_SBG_UaW
  CASE('X_LES_RES_ddz_Rt_SBG_W2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddz_Rt_SBG_W2,NDIMS)
    X_LES_RES_ddz_Rt_SBG_W2=>TLES%X_LES_RES_ddz_Rt_SBG_W2
  CASE('X_LES_RES_ddxa_Thl_SBG_UaRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Thl_SBG_UaRt,NDIMS)
    X_LES_RES_ddxa_Thl_SBG_UaRt=>TLES%X_LES_RES_ddxa_Thl_SBG_UaRt
  CASE('X_LES_RES_ddxa_Rt_SBG_UaThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Rt_SBG_UaThl,NDIMS)
    X_LES_RES_ddxa_Rt_SBG_UaThl=>TLES%X_LES_RES_ddxa_Rt_SBG_UaThl
  CASE('X_LES_RES_ddxa_Thl_SBG_UaThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Thl_SBG_UaThl,NDIMS)
    X_LES_RES_ddxa_Thl_SBG_UaThl=>TLES%X_LES_RES_ddxa_Thl_SBG_UaThl
  CASE('X_LES_RES_ddxa_Rt_SBG_UaRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Rt_SBG_UaRt,NDIMS)
    X_LES_RES_ddxa_Rt_SBG_UaRt=>TLES%X_LES_RES_ddxa_Rt_SBG_UaRt
  CASE('X_LES_RES_ddxa_W_SBG_UaSv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_W_SBG_UaSv,NDIMS)
    X_LES_RES_ddxa_W_SBG_UaSv=>TLES%X_LES_RES_ddxa_W_SBG_UaSv
  CASE('X_LES_RES_ddxa_Sv_SBG_UaW')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Sv_SBG_UaW,NDIMS)
    X_LES_RES_ddxa_Sv_SBG_UaW=>TLES%X_LES_RES_ddxa_Sv_SBG_UaW
  CASE('X_LES_RES_ddz_Sv_SBG_W2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddz_Sv_SBG_W2,NDIMS)
    X_LES_RES_ddz_Sv_SBG_W2=>TLES%X_LES_RES_ddz_Sv_SBG_W2
  CASE('X_LES_RES_ddxa_Sv_SBG_UaSv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_RES_ddxa_Sv_SBG_UaSv,NDIMS)
    X_LES_RES_ddxa_Sv_SBG_UaSv=>TLES%X_LES_RES_ddxa_Sv_SBG_UaSv
  CASE('X_LES_SUBGRID_U2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_U2,NDIMS)
    X_LES_SUBGRID_U2=>TLES%X_LES_SUBGRID_U2
  CASE('X_LES_SUBGRID_V2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_V2,NDIMS)
    X_LES_SUBGRID_V2=>TLES%X_LES_SUBGRID_V2
  CASE('X_LES_SUBGRID_W2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_W2,NDIMS)
    X_LES_SUBGRID_W2=>TLES%X_LES_SUBGRID_W2
  CASE('X_LES_SUBGRID_Thl2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Thl2,NDIMS)
    X_LES_SUBGRID_Thl2=>TLES%X_LES_SUBGRID_Thl2
  CASE('X_LES_SUBGRID_Rt2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Rt2,NDIMS)
    X_LES_SUBGRID_Rt2=>TLES%X_LES_SUBGRID_Rt2
  CASE('X_LES_SUBGRID_Rc2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Rc2,NDIMS)
    X_LES_SUBGRID_Rc2=>TLES%X_LES_SUBGRID_Rc2
  CASE('X_LES_SUBGRID_Ri2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Ri2,NDIMS)
    X_LES_SUBGRID_Ri2=>TLES%X_LES_SUBGRID_Ri2
  CASE('X_LES_SUBGRID_ThlRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_ThlRt,NDIMS)
    X_LES_SUBGRID_ThlRt=>TLES%X_LES_SUBGRID_ThlRt
  CASE('X_LES_SUBGRID_Sv2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Sv2,NDIMS)
    X_LES_SUBGRID_Sv2=>TLES%X_LES_SUBGRID_Sv2
  CASE('X_LES_SUBGRID_UV')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_UV,NDIMS)
    X_LES_SUBGRID_UV=>TLES%X_LES_SUBGRID_UV
  CASE('X_LES_SUBGRID_WU')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WU,NDIMS)
    X_LES_SUBGRID_WU=>TLES%X_LES_SUBGRID_WU
  CASE('X_LES_SUBGRID_WV')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WV,NDIMS)
    X_LES_SUBGRID_WV=>TLES%X_LES_SUBGRID_WV
  CASE('X_LES_SUBGRID_UThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_UThl,NDIMS)
    X_LES_SUBGRID_UThl=>TLES%X_LES_SUBGRID_UThl
  CASE('X_LES_SUBGRID_VThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_VThl,NDIMS)
    X_LES_SUBGRID_VThl=>TLES%X_LES_SUBGRID_VThl
  CASE('X_LES_SUBGRID_WThl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WThl,NDIMS)
    X_LES_SUBGRID_WThl=>TLES%X_LES_SUBGRID_WThl
  CASE('X_LES_SUBGRID_URt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_URt,NDIMS)
    X_LES_SUBGRID_URt=>TLES%X_LES_SUBGRID_URt
  CASE('X_LES_SUBGRID_VRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_VRt,NDIMS)
    X_LES_SUBGRID_VRt=>TLES%X_LES_SUBGRID_VRt
  CASE('X_LES_SUBGRID_WRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WRt,NDIMS)
    X_LES_SUBGRID_WRt=>TLES%X_LES_SUBGRID_WRt
  CASE('X_LES_SUBGRID_URc')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_URc,NDIMS)
    X_LES_SUBGRID_URc=>TLES%X_LES_SUBGRID_URc
  CASE('X_LES_SUBGRID_VRc')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_VRc,NDIMS)
    X_LES_SUBGRID_VRc=>TLES%X_LES_SUBGRID_VRc
  CASE('X_LES_SUBGRID_WRc')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WRc,NDIMS)
    X_LES_SUBGRID_WRc=>TLES%X_LES_SUBGRID_WRc
  CASE('X_LES_SUBGRID_USv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_USv,NDIMS)
    X_LES_SUBGRID_USv=>TLES%X_LES_SUBGRID_USv
  CASE('X_LES_SUBGRID_VSv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_VSv,NDIMS)
    X_LES_SUBGRID_VSv=>TLES%X_LES_SUBGRID_VSv
  CASE('X_LES_SUBGRID_WSv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WSv,NDIMS)
    X_LES_SUBGRID_WSv=>TLES%X_LES_SUBGRID_WSv
  CASE('X_LES_SUBGRID_UTke')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_UTke,NDIMS)
    X_LES_SUBGRID_UTke=>TLES%X_LES_SUBGRID_UTke
  CASE('X_LES_SUBGRID_VTke')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_VTke,NDIMS)
    X_LES_SUBGRID_VTke=>TLES%X_LES_SUBGRID_VTke
  CASE('X_LES_SUBGRID_WTke')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WTke,NDIMS)
    X_LES_SUBGRID_WTke=>TLES%X_LES_SUBGRID_WTke
  CASE('X_LES_SUBGRID_ddz_WTke')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_ddz_WTke,NDIMS)
    X_LES_SUBGRID_ddz_WTke=>TLES%X_LES_SUBGRID_ddz_WTke
  CASE('X_LES_SUBGRID_WThv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WThv,NDIMS)
    X_LES_SUBGRID_WThv=>TLES%X_LES_SUBGRID_WThv
  CASE('X_LES_SUBGRID_ThlThv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_ThlThv,NDIMS)
    X_LES_SUBGRID_ThlThv=>TLES%X_LES_SUBGRID_ThlThv
  CASE('X_LES_SUBGRID_RtThv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RtThv,NDIMS)
    X_LES_SUBGRID_RtThv=>TLES%X_LES_SUBGRID_RtThv
  CASE('X_LES_SUBGRID_SvThv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_SvThv,NDIMS)
    X_LES_SUBGRID_SvThv=>TLES%X_LES_SUBGRID_SvThv
  CASE('X_LES_SUBGRID_W2Thl')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_W2Thl,NDIMS)
    X_LES_SUBGRID_W2Thl=>TLES%X_LES_SUBGRID_W2Thl
  CASE('X_LES_SUBGRID_W2Rt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_W2Rt,NDIMS)
    X_LES_SUBGRID_W2Rt=>TLES%X_LES_SUBGRID_W2Rt
  CASE('X_LES_SUBGRID_W2Sv')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_W2Sv,NDIMS)
    X_LES_SUBGRID_W2Sv=>TLES%X_LES_SUBGRID_W2Sv
  CASE('X_LES_SUBGRID_WThlRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WThlRt,NDIMS)
    X_LES_SUBGRID_WThlRt=>TLES%X_LES_SUBGRID_WThlRt
  CASE('X_LES_SUBGRID_WThl2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WThl2,NDIMS)
    X_LES_SUBGRID_WThl2=>TLES%X_LES_SUBGRID_WThl2
  CASE('X_LES_SUBGRID_WRt2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WRt2,NDIMS)
    X_LES_SUBGRID_WRt2=>TLES%X_LES_SUBGRID_WRt2
  CASE('X_LES_SUBGRID_WSv2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WSv2,NDIMS)
    X_LES_SUBGRID_WSv2=>TLES%X_LES_SUBGRID_WSv2
  CASE('X_LES_SUBGRID_DISS_Tke')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DISS_Tke,NDIMS)
    X_LES_SUBGRID_DISS_Tke=>TLES%X_LES_SUBGRID_DISS_Tke
  CASE('X_LES_SUBGRID_DISS_Thl2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DISS_Thl2,NDIMS)
    X_LES_SUBGRID_DISS_Thl2=>TLES%X_LES_SUBGRID_DISS_Thl2
  CASE('X_LES_SUBGRID_DISS_Rt2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DISS_Rt2,NDIMS)
    X_LES_SUBGRID_DISS_Rt2=>TLES%X_LES_SUBGRID_DISS_Rt2
  CASE('X_LES_SUBGRID_DISS_ThlRt')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DISS_ThlRt,NDIMS)
    X_LES_SUBGRID_DISS_ThlRt=>TLES%X_LES_SUBGRID_DISS_ThlRt
  CASE('X_LES_SUBGRID_DISS_Sv2')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DISS_Sv2,NDIMS)
    X_LES_SUBGRID_DISS_Sv2=>TLES%X_LES_SUBGRID_DISS_Sv2
  CASE('X_LES_SUBGRID_WP')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WP,NDIMS)
    X_LES_SUBGRID_WP=>TLES%X_LES_SUBGRID_WP
  CASE('X_LES_SUBGRID_ThlPz')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_ThlPz,NDIMS)
    X_LES_SUBGRID_ThlPz=>TLES%X_LES_SUBGRID_ThlPz
  CASE('X_LES_SUBGRID_RtPz')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RtPz,NDIMS)
    X_LES_SUBGRID_RtPz=>TLES%X_LES_SUBGRID_RtPz
  CASE('X_LES_SUBGRID_SvPz')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_SvPz,NDIMS)
    X_LES_SUBGRID_SvPz=>TLES%X_LES_SUBGRID_SvPz
  CASE('X_LES_SUBGRID_PHI3')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_PHI3,NDIMS)
    X_LES_SUBGRID_PHI3=>TLES%X_LES_SUBGRID_PHI3
  CASE('X_LES_SUBGRID_PSI3')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_PSI3,NDIMS)
    X_LES_SUBGRID_PSI3=>TLES%X_LES_SUBGRID_PSI3
  CASE('X_LES_SUBGRID_LMix')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_LMix,NDIMS)
    X_LES_SUBGRID_LMix=>TLES%X_LES_SUBGRID_LMix
  CASE('X_LES_SUBGRID_LDiss')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_LDiss,NDIMS)
    X_LES_SUBGRID_LDiss=>TLES%X_LES_SUBGRID_LDiss
  CASE('X_LES_SUBGRID_Km')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Km,NDIMS)
    X_LES_SUBGRID_Km=>TLES%X_LES_SUBGRID_Km
  CASE('X_LES_SUBGRID_Kh')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_Kh,NDIMS)
    X_LES_SUBGRID_Kh=>TLES%X_LES_SUBGRID_Kh
  CASE('X_LES_SUBGRID_THLUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_THLUP_MF,NDIMS)
    X_LES_SUBGRID_THLUP_MF=>TLES%X_LES_SUBGRID_THLUP_MF
  CASE('X_LES_SUBGRID_RTUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RTUP_MF,NDIMS)
    X_LES_SUBGRID_RTUP_MF=>TLES%X_LES_SUBGRID_RTUP_MF
  CASE('X_LES_SUBGRID_RVUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RVUP_MF,NDIMS)
    X_LES_SUBGRID_RVUP_MF=>TLES%X_LES_SUBGRID_RVUP_MF
  CASE('X_LES_SUBGRID_RCUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RCUP_MF,NDIMS)
    X_LES_SUBGRID_RCUP_MF=>TLES%X_LES_SUBGRID_RCUP_MF
  CASE('X_LES_SUBGRID_RIUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_RIUP_MF,NDIMS)
    X_LES_SUBGRID_RIUP_MF=>TLES%X_LES_SUBGRID_RIUP_MF
  CASE('X_LES_SUBGRID_WUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WUP_MF,NDIMS)
    X_LES_SUBGRID_WUP_MF=>TLES%X_LES_SUBGRID_WUP_MF
  CASE('X_LES_SUBGRID_MASSFLUX')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_MASSFLUX,NDIMS)
    X_LES_SUBGRID_MASSFLUX=>TLES%X_LES_SUBGRID_MASSFLUX
  CASE('X_LES_SUBGRID_DETR')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_DETR,NDIMS)
    X_LES_SUBGRID_DETR=>TLES%X_LES_SUBGRID_DETR
  CASE('X_LES_SUBGRID_ENTR')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_ENTR,NDIMS)
    X_LES_SUBGRID_ENTR=>TLES%X_LES_SUBGRID_ENTR
  CASE('X_LES_SUBGRID_FRACUP')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_FRACUP,NDIMS)
    X_LES_SUBGRID_FRACUP=>TLES%X_LES_SUBGRID_FRACUP
  CASE('X_LES_SUBGRID_THVUP_MF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_THVUP_MF,NDIMS)
    X_LES_SUBGRID_THVUP_MF=>TLES%X_LES_SUBGRID_THVUP_MF
  CASE('X_LES_SUBGRID_WTHLMF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WTHLMF,NDIMS)
    X_LES_SUBGRID_WTHLMF=>TLES%X_LES_SUBGRID_WTHLMF
  CASE('X_LES_SUBGRID_WRTMF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WRTMF,NDIMS)
    X_LES_SUBGRID_WRTMF=>TLES%X_LES_SUBGRID_WRTMF
  CASE('X_LES_SUBGRID_WTHVMF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WTHVMF,NDIMS)
    X_LES_SUBGRID_WTHVMF=>TLES%X_LES_SUBGRID_WTHVMF
  CASE('X_LES_SUBGRID_WUMF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WUMF,NDIMS)
    X_LES_SUBGRID_WUMF=>TLES%X_LES_SUBGRID_WUMF
  CASE('X_LES_SUBGRID_WVMF')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SUBGRID_WVMF,NDIMS)
    X_LES_SUBGRID_WVMF=>TLES%X_LES_SUBGRID_WVMF
  CASE('X_LES_USTAR')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_USTAR,NDIMS)
    X_LES_USTAR=>TLES%X_LES_USTAR
  CASE('X_LES_UW0')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_UW0,NDIMS)
    X_LES_UW0=>TLES%X_LES_UW0
  CASE('X_LES_VW0')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_VW0,NDIMS)
    X_LES_VW0=>TLES%X_LES_VW0
  CASE('X_LES_Q0')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_Q0,NDIMS)
    X_LES_Q0=>TLES%X_LES_Q0
  CASE('X_LES_E0')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_E0,NDIMS)
    X_LES_E0=>TLES%X_LES_E0
  CASE('X_LES_SV0')
    CALL LES_ALLOCATE_DIM(TLES%X_LES_SV0,NDIMS)
    X_LES_SV0=>TLES%X_LES_SV0
  CASE('XLES_PDF_RV')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RV,NDIMS)
    XLES_PDF_RV=>TLES%XLES_PDF_RV
  CASE('XLES_PDF_TH')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_TH,NDIMS)
    XLES_PDF_TH=>TLES%XLES_PDF_TH
  CASE('XLES_PDF_W')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_W,NDIMS)
    XLES_PDF_W=>TLES%XLES_PDF_W
  CASE('XLES_PDF_THV')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_THV,NDIMS)
    XLES_PDF_THV=>TLES%XLES_PDF_THV
  CASE('XLES_PDF_RC')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RC,NDIMS)
    XLES_PDF_RC=>TLES%XLES_PDF_RC
  CASE('XLES_PDF_RR')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RR,NDIMS)
    XLES_PDF_RR=>TLES%XLES_PDF_RR
  CASE('XLES_PDF_RI')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RI,NDIMS)
    XLES_PDF_RI=>TLES%XLES_PDF_RI
  CASE('XLES_PDF_RS')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RS,NDIMS)
    XLES_PDF_RS=>TLES%XLES_PDF_RS
  CASE('XLES_PDF_RG')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RG,NDIMS)
    XLES_PDF_RG=>TLES%XLES_PDF_RG
  CASE('XLES_PDF_RT')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_RT,NDIMS)
    XLES_PDF_RT=>TLES%XLES_PDF_RT
  CASE('XLES_PDF_THL')
    CALL LES_ALLOCATE_DIM(TLES%XLES_PDF_THL,NDIMS)
    XLES_PDF_THL=>TLES%XLES_PDF_THL
  END SELECT
  !
END SUBROUTINE LES_ALLOCATE
!
SUBROUTINE LES_DEALLOCATE(HNAME)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
!
  SELECT CASE(HNAME)
  CASE('LLES_CURRENT_CART_MASK')
    LLES_CURRENT_CART_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_CART_MASK)
  CASE('LLES_CURRENT_NEB_MASK')
    LLES_CURRENT_NEB_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_NEB_MASK)
  CASE('LLES_CURRENT_CORE_MASK')
    LLES_CURRENT_CORE_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_CORE_MASK)
  CASE('LLES_CURRENT_MY_MASKS')
    LLES_CURRENT_MY_MASKS=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_MY_MASKS)
  CASE('LLES_CURRENT_CS1_MASK')
    LLES_CURRENT_CS1_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_CS1_MASK)
  CASE('LLES_CURRENT_CS2_MASK')
    LLES_CURRENT_CS2_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_CS2_MASK)
  CASE('LLES_CURRENT_CS3_MASK')
    LLES_CURRENT_CS3_MASK=>NULL()
    DEALLOCATE(TLES%LLES_CURRENT_CS3_MASK)
  CASE('XLES_CURRENT_Z')
    XLES_CURRENT_Z=>NULL()
    DEALLOCATE(TLES%XLES_CURRENT_Z)
  CASE('NKLIN_CURRENT_LES')
    NKLIN_CURRENT_LES=>NULL()
    DEALLOCATE(TLES%NKLIN_CURRENT_LES)
  CASE('XCOEFLIN_CURRENT_LES')
    XCOEFLIN_CURRENT_LES=>NULL()
    DEALLOCATE(TLES%XCOEFLIN_CURRENT_LES)
  CASE('NKLIN_CURRENT_SPEC')
    NKLIN_CURRENT_SPEC=>NULL()
    DEALLOCATE(TLES%NKLIN_CURRENT_SPEC)
  CASE('XCOEFLIN_CURRENT_SPEC')
    XCOEFLIN_CURRENT_SPEC=>NULL()
    DEALLOCATE(TLES%XCOEFLIN_CURRENT_SPEC)
  CASE('XLES_NORM_M')
    XLES_NORM_M=>NULL()
    DEALLOCATE(TLES%XLES_NORM_M)
  CASE('XLES_NORM_K')
    XLES_NORM_K=>NULL()
    DEALLOCATE(TLES%XLES_NORM_K)
  CASE('XLES_NORM_S')
    XLES_NORM_S=>NULL()
    DEALLOCATE(TLES%XLES_NORM_S)
  CASE('XLES_NORM_RHO')
    XLES_NORM_RHO=>NULL()
    DEALLOCATE(TLES%XLES_NORM_RHO)
  CASE('XLES_NORM_RV')
    XLES_NORM_RV=>NULL()
    DEALLOCATE(TLES%XLES_NORM_RV)
  CASE('XLES_NORM_SV')
    XLES_NORM_SV=>NULL()
    DEALLOCATE(TLES%XLES_NORM_SV)
  CASE('XLES_NORM_P')
    XLES_NORM_P=>NULL()
    DEALLOCATE(TLES%XLES_NORM_P)
  CASE('X_LES_RES_W_SBG_WThl')
    X_LES_RES_W_SBG_WThl=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_WThl)
  CASE('X_LES_RES_W_SBG_WRt')
    X_LES_RES_W_SBG_WRt=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_WRt)
  CASE('X_LES_RES_W_SBG_Thl2')
    X_LES_RES_W_SBG_Thl2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_Thl2)
  CASE('X_LES_RES_W_SBG_Rt2')
    X_LES_RES_W_SBG_Rt2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_Rt2)
  CASE('X_LES_RES_W_SBG_ThlRt')
    X_LES_RES_W_SBG_ThlRt=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_ThlRt)
  CASE('X_LES_RES_W_SBG_WSv')
    X_LES_RES_W_SBG_WSv=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_WSv)
  CASE('X_LES_RES_W_SBG_Sv2')
    X_LES_RES_W_SBG_Sv2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_W_SBG_Sv2)
  CASE('XLES_SUBGRID_RCSIGS')
    XLES_SUBGRID_RCSIGS=>NULL()
    DEALLOCATE(TLES%XLES_SUBGRID_RCSIGS)
  CASE('XLES_SUBGRID_RCSIGC')
    XLES_SUBGRID_RCSIGC=>NULL()
    DEALLOCATE(TLES%XLES_SUBGRID_RCSIGC)
  CASE('X_LES_RES_ddxa_U_SBG_UaU')
    X_LES_RES_ddxa_U_SBG_UaU=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_U_SBG_UaU)
  CASE('X_LES_RES_ddxa_V_SBG_UaV')
    X_LES_RES_ddxa_V_SBG_UaV=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_V_SBG_UaV)
  CASE('X_LES_RES_ddxa_W_SBG_UaW')
    X_LES_RES_ddxa_W_SBG_UaW=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_W_SBG_UaW)
  CASE('X_LES_RES_ddxa_W_SBG_UaThl')
    X_LES_RES_ddxa_W_SBG_UaThl=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_W_SBG_UaThl)
  CASE('X_LES_RES_ddxa_Thl_SBG_UaW')
    X_LES_RES_ddxa_Thl_SBG_UaW=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Thl_SBG_UaW)
  CASE('X_LES_RES_ddz_Thl_SBG_W2')
    X_LES_RES_ddz_Thl_SBG_W2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddz_Thl_SBG_W2)
  CASE('X_LES_RES_ddxa_W_SBG_UaRt')
    X_LES_RES_ddxa_W_SBG_UaRt=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_W_SBG_UaRt)
  CASE('X_LES_RES_ddxa_Rt_SBG_UaW')
    X_LES_RES_ddxa_Rt_SBG_UaW=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Rt_SBG_UaW)
  CASE('X_LES_RES_ddz_Rt_SBG_W2')
    X_LES_RES_ddz_Rt_SBG_W2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddz_Rt_SBG_W2)
  CASE('X_LES_RES_ddxa_Thl_SBG_UaRt')
    X_LES_RES_ddxa_Thl_SBG_UaRt=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Thl_SBG_UaRt)
  CASE('X_LES_RES_ddxa_Rt_SBG_UaThl')
    X_LES_RES_ddxa_Rt_SBG_UaThl=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Rt_SBG_UaThl)
  CASE('X_LES_RES_ddxa_Thl_SBG_UaThl')
    X_LES_RES_ddxa_Thl_SBG_UaThl=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Thl_SBG_UaThl)
  CASE('X_LES_RES_ddxa_Rt_SBG_UaRt')
    X_LES_RES_ddxa_Rt_SBG_UaRt=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Rt_SBG_UaRt)
  CASE('X_LES_RES_ddxa_W_SBG_UaSv')
    X_LES_RES_ddxa_W_SBG_UaSv=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_W_SBG_UaSv)
  CASE('X_LES_RES_ddxa_Sv_SBG_UaW')
    X_LES_RES_ddxa_Sv_SBG_UaW=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Sv_SBG_UaW)
  CASE('X_LES_RES_ddz_Sv_SBG_W2')
    X_LES_RES_ddz_Sv_SBG_W2=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddz_Sv_SBG_W2)
  CASE('X_LES_RES_ddxa_Sv_SBG_UaSv')
    X_LES_RES_ddxa_Sv_SBG_UaSv=>NULL()
    DEALLOCATE(TLES%X_LES_RES_ddxa_Sv_SBG_UaSv)
  CASE('X_LES_SUBGRID_U2')
    X_LES_SUBGRID_U2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_U2)
  CASE('X_LES_SUBGRID_V2')
    X_LES_SUBGRID_V2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_V2)
  CASE('X_LES_SUBGRID_W2')
    X_LES_SUBGRID_W2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_W2)
  CASE('X_LES_SUBGRID_Thl2')
    X_LES_SUBGRID_Thl2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Thl2)
  CASE('X_LES_SUBGRID_Rt2')
    X_LES_SUBGRID_Rt2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Rt2)
  CASE('X_LES_SUBGRID_Rc2')
    X_LES_SUBGRID_Rc2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Rc2)
  CASE('X_LES_SUBGRID_Ri2')
    X_LES_SUBGRID_Ri2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Ri2)
  CASE('X_LES_SUBGRID_ThlRt')
    X_LES_SUBGRID_ThlRt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_ThlRt)
  CASE('X_LES_SUBGRID_Sv2')
    X_LES_SUBGRID_Sv2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Sv2)
  CASE('X_LES_SUBGRID_UV')
    X_LES_SUBGRID_UV=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_UV)
  CASE('X_LES_SUBGRID_WU')
    X_LES_SUBGRID_WU=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WU)
  CASE('X_LES_SUBGRID_WV')
    X_LES_SUBGRID_WV=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WV)
  CASE('X_LES_SUBGRID_UThl')
    X_LES_SUBGRID_UThl=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_UThl)
  CASE('X_LES_SUBGRID_VThl')
    X_LES_SUBGRID_VThl=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_VThl)
  CASE('X_LES_SUBGRID_WThl')
    X_LES_SUBGRID_WThl=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WThl)
  CASE('X_LES_SUBGRID_URt')
    X_LES_SUBGRID_URt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_URt)
  CASE('X_LES_SUBGRID_VRt')
    X_LES_SUBGRID_VRt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_VRt)
  CASE('X_LES_SUBGRID_WRt')
    X_LES_SUBGRID_WRt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WRt)
  CASE('X_LES_SUBGRID_URc')
    X_LES_SUBGRID_URc=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_URc)
  CASE('X_LES_SUBGRID_VRc')
    X_LES_SUBGRID_VRc=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_VRc)
  CASE('X_LES_SUBGRID_WRc')
    X_LES_SUBGRID_WRc=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WRc)
  CASE('X_LES_SUBGRID_USv')
    X_LES_SUBGRID_USv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_USv)
  CASE('X_LES_SUBGRID_VSv')
    X_LES_SUBGRID_VSv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_VSv)
  CASE('X_LES_SUBGRID_WSv')
    X_LES_SUBGRID_WSv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WSv)
  CASE('X_LES_SUBGRID_UTke')
    X_LES_SUBGRID_UTke=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_UTke)
  CASE('X_LES_SUBGRID_VTke')
    X_LES_SUBGRID_VTke=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_VTke)
  CASE('X_LES_SUBGRID_WTke')
    X_LES_SUBGRID_WTke=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WTke)
  CASE('X_LES_SUBGRID_ddz_WTke')
    X_LES_SUBGRID_ddz_WTke=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_ddz_WTke)
  CASE('X_LES_SUBGRID_WThv')
    X_LES_SUBGRID_WThv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WThv)
  CASE('X_LES_SUBGRID_ThlThv')
    X_LES_SUBGRID_ThlThv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_ThlThv)
  CASE('X_LES_SUBGRID_RtThv')
    X_LES_SUBGRID_RtThv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RtThv)
  CASE('X_LES_SUBGRID_SvThv')
    X_LES_SUBGRID_SvThv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_SvThv)
  CASE('X_LES_SUBGRID_W2Thl')
    X_LES_SUBGRID_W2Thl=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_W2Thl)
  CASE('X_LES_SUBGRID_W2Rt')
    X_LES_SUBGRID_W2Rt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_W2Rt)
  CASE('X_LES_SUBGRID_W2Sv')
    X_LES_SUBGRID_W2Sv=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_W2Sv)
  CASE('X_LES_SUBGRID_WThlRt')
    X_LES_SUBGRID_WThlRt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WThlRt)
  CASE('X_LES_SUBGRID_WThl2')
    X_LES_SUBGRID_WThl2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WThl2)
  CASE('X_LES_SUBGRID_WRt2')
    X_LES_SUBGRID_WRt2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WRt2)
  CASE('X_LES_SUBGRID_WSv2')
    X_LES_SUBGRID_WSv2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WSv2)
  CASE('X_LES_SUBGRID_DISS_Tke')
    X_LES_SUBGRID_DISS_Tke=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DISS_Tke)
  CASE('X_LES_SUBGRID_DISS_Thl2')
    X_LES_SUBGRID_DISS_Thl2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DISS_Thl2)
  CASE('X_LES_SUBGRID_DISS_Rt2')
    X_LES_SUBGRID_DISS_Rt2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DISS_Rt2)
  CASE('X_LES_SUBGRID_DISS_ThlRt')
    X_LES_SUBGRID_DISS_ThlRt=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DISS_ThlRt)
  CASE('X_LES_SUBGRID_DISS_Sv2')
    X_LES_SUBGRID_DISS_Sv2=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DISS_Sv2)
  CASE('X_LES_SUBGRID_WP')
    X_LES_SUBGRID_WP=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WP)
  CASE('X_LES_SUBGRID_ThlPz')
    X_LES_SUBGRID_ThlPz=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_ThlPz)
  CASE('X_LES_SUBGRID_RtPz')
    X_LES_SUBGRID_RtPz=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RtPz)
  CASE('X_LES_SUBGRID_SvPz')
    X_LES_SUBGRID_SvPz=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_SvPz)
  CASE('X_LES_SUBGRID_PHI3')
    X_LES_SUBGRID_PHI3=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_PHI3)
  CASE('X_LES_SUBGRID_PSI3')
    X_LES_SUBGRID_PSI3=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_PSI3)
  CASE('X_LES_SUBGRID_LMix')
    X_LES_SUBGRID_LMix=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_LMix)
  CASE('X_LES_SUBGRID_LDiss')
    X_LES_SUBGRID_LDiss=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_LDiss)
  CASE('X_LES_SUBGRID_Km')
    X_LES_SUBGRID_Km=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Km)
  CASE('X_LES_SUBGRID_Kh')
    X_LES_SUBGRID_Kh=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_Kh)
  CASE('X_LES_SUBGRID_THLUP_MF')
    X_LES_SUBGRID_THLUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_THLUP_MF)
  CASE('X_LES_SUBGRID_RTUP_MF')
    X_LES_SUBGRID_RTUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RTUP_MF)
  CASE('X_LES_SUBGRID_RVUP_MF')
    X_LES_SUBGRID_RVUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RVUP_MF)
  CASE('X_LES_SUBGRID_RCUP_MF')
    X_LES_SUBGRID_RCUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RCUP_MF)
  CASE('X_LES_SUBGRID_RIUP_MF')
    X_LES_SUBGRID_RIUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_RIUP_MF)
  CASE('X_LES_SUBGRID_WUP_MF')
    X_LES_SUBGRID_WUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WUP_MF)
  CASE('X_LES_SUBGRID_MASSFLUX')
    X_LES_SUBGRID_MASSFLUX=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_MASSFLUX)
  CASE('X_LES_SUBGRID_DETR')
    X_LES_SUBGRID_DETR=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_DETR)
  CASE('X_LES_SUBGRID_ENTR')
    X_LES_SUBGRID_ENTR=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_ENTR)
  CASE('X_LES_SUBGRID_FRACUP')
    X_LES_SUBGRID_FRACUP=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_FRACUP)
  CASE('X_LES_SUBGRID_THVUP_MF')
    X_LES_SUBGRID_THVUP_MF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_THVUP_MF)
  CASE('X_LES_SUBGRID_WTHLMF')
    X_LES_SUBGRID_WTHLMF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WTHLMF)
  CASE('X_LES_SUBGRID_WRTMF')
    X_LES_SUBGRID_WRTMF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WRTMF)
  CASE('X_LES_SUBGRID_WTHVMF')
    X_LES_SUBGRID_WTHVMF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WTHVMF)
  CASE('X_LES_SUBGRID_WUMF')
    X_LES_SUBGRID_WUMF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WUMF)
  CASE('X_LES_SUBGRID_WVMF')
    X_LES_SUBGRID_WVMF=>NULL()
    DEALLOCATE(TLES%X_LES_SUBGRID_WVMF)
  CASE('X_LES_USTAR')
    X_LES_USTAR=>NULL()
    DEALLOCATE(TLES%X_LES_USTAR)
  CASE('X_LES_UW0')
    X_LES_UW0=>NULL()
    DEALLOCATE(TLES%X_LES_UW0)
  CASE('X_LES_VW0')
    X_LES_VW0=>NULL()
    DEALLOCATE(TLES%X_LES_VW0)
  CASE('X_LES_Q0')
    X_LES_Q0=>NULL()
    DEALLOCATE(TLES%X_LES_Q0)
  CASE('X_LES_E0')
    X_LES_E0=>NULL()
    DEALLOCATE(TLES%X_LES_E0)
  CASE('X_LES_SV0')
    X_LES_SV0=>NULL()
    DEALLOCATE(TLES%X_LES_SV0)
  CASE('XLES_PDF_RV')
    XLES_PDF_RV=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RV)
  CASE('XLES_PDF_TH')
    XLES_PDF_TH=>NULL()
    DEALLOCATE(TLES%XLES_PDF_TH)
  CASE('XLES_PDF_W')
    XLES_PDF_W=>NULL()
    DEALLOCATE(TLES%XLES_PDF_W)
  CASE('XLES_PDF_THV')
    XLES_PDF_THV=>NULL()
    DEALLOCATE(TLES%XLES_PDF_THV)
  CASE('XLES_PDF_RC')
    XLES_PDF_RC=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RC)
  CASE('XLES_PDF_RR')
    XLES_PDF_RR=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RR)
  CASE('XLES_PDF_RI')
    XLES_PDF_RI=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RI)
  CASE('XLES_PDF_RS')
    XLES_PDF_RS=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RS)
  CASE('XLES_PDF_RG')
    XLES_PDF_RG=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RG)
  CASE('XLES_PDF_RT')
    XLES_PDF_RT=>NULL()
    DEALLOCATE(TLES%XLES_PDF_RT)
  CASE('XLES_PDF_THL')
    XLES_PDF_THL=>NULL()
    DEALLOCATE(TLES%XLES_PDF_THL)
  END SELECT
END SUBROUTINE LES_DEALLOCATE
!!
!SUBROUTINE LES_INI_TIMESTEP_DEALLOCATE()
!  IMPLICIT NONE
!  XCOEFLIN_CURRENT_SPEC=>NULL()
!  DEALLOCATE(TLES%XCOEFLIN_CURRENT_SPEC)
!END SUBROUTINE LES_INI_TIMESTEP_DEALLOCATE
!
SUBROUTINE LES_ALLOCATE_1DIMX(PVAR,KDIM)
  IMPLICIT NONE
  REAL, DIMENSION(:),ALLOCATABLE,  INTENT(OUT) :: PVAR
  INTEGER, DIMENSION(1), INTENT(IN) :: KDIM
  ALLOCATE(PVAR(KDIM(1)))
END SUBROUTINE LES_ALLOCATE_1DIMX
!
SUBROUTINE LES_ALLOCATE_2DIMX(PVAR,KDIM)
  IMPLICIT NONE
  REAL, DIMENSION(:,:),ALLOCATABLE, INTENT(OUT) :: PVAR
  INTEGER, DIMENSION(2), INTENT(IN) :: KDIM
  ALLOCATE(PVAR(KDIM(1),KDIM(2)))
END SUBROUTINE LES_ALLOCATE_2DIMX
!
SUBROUTINE LES_ALLOCATE_3DIMX(PVAR,KDIM)
  IMPLICIT NONE
  REAL, DIMENSION(:,:,:),ALLOCATABLE,  INTENT(OUT) :: PVAR
  INTEGER, DIMENSION(3), INTENT(IN) :: KDIM
  ALLOCATE(PVAR(KDIM(1),KDIM(2),KDIM(3)))
END SUBROUTINE LES_ALLOCATE_3DIMX
!
SUBROUTINE LES_ALLOCATE_4DIMX(PVAR,KDIM)
  IMPLICIT NONE
  REAL, DIMENSION(:,:,:,:),ALLOCATABLE,  INTENT(OUT) :: PVAR
  INTEGER, DIMENSION(4), INTENT(IN) :: KDIM
  ALLOCATE(PVAR(KDIM(1),KDIM(2),KDIM(3),KDIM(4)))
END SUBROUTINE LES_ALLOCATE_4DIMX
!
SUBROUTINE LES_ALLOCATE_1DIMI(KVAR,KDIM)
  IMPLICIT NONE
  INTEGER, DIMENSION(:),ALLOCATABLE,  INTENT(OUT) :: KVAR
  INTEGER, DIMENSION(1), INTENT(IN) :: KDIM
  ALLOCATE(KVAR(KDIM(1)))
END SUBROUTINE LES_ALLOCATE_1DIMI
!
SUBROUTINE LES_ALLOCATE_3DIMI(KVAR,KDIM)
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:,:),ALLOCATABLE,  INTENT(OUT) :: KVAR
  INTEGER, DIMENSION(3), INTENT(IN) :: KDIM
  ALLOCATE(KVAR(KDIM(1),KDIM(2),KDIM(3)))
END SUBROUTINE LES_ALLOCATE_3DIMI
!
SUBROUTINE LES_ALLOCATE_3DIML(OVAR,KDIM)
  IMPLICIT NONE
  LOGICAL, DIMENSION(:,:,:),ALLOCATABLE,  INTENT(OUT) :: OVAR
  INTEGER, DIMENSION(3), INTENT(IN) :: KDIM
  ALLOCATE(OVAR(KDIM(1),KDIM(2),KDIM(3)))
END SUBROUTINE LES_ALLOCATE_3DIML
!
SUBROUTINE LES_ALLOCATE_4DIML(OVAR,KDIM)
  IMPLICIT NONE
  LOGICAL, DIMENSION(:,:,:,:),ALLOCATABLE,  INTENT(OUT) :: OVAR
  INTEGER, DIMENSION(4), INTENT(IN) :: KDIM
  ALLOCATE(OVAR(KDIM(1),KDIM(2),KDIM(3),KDIM(4)))
END SUBROUTINE LES_ALLOCATE_4DIML
!
SUBROUTINE LES_ALLOCATE_2DIMC(HVAR,KDIM)
  IMPLICIT NONE
  LOGICAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: HVAR
  INTEGER, DIMENSION(2), INTENT(IN) :: KDIM
  ALLOCATE(HVAR(KDIM(1),KDIM(2)))
END SUBROUTINE LES_ALLOCATE_2DIMC
!
END MODULE MODD_LES
