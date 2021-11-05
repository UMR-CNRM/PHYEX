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
!!	   J. Cuxart   *INM and Meteo France*
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
INTEGER, DIMENSION( 10) :: NLES_TEMP_SERIE_I   ! I, J and Z point
INTEGER, DIMENSION( 10) :: NLES_TEMP_SERIE_J   ! localizations to
INTEGER, DIMENSION( 10) :: NLES_TEMP_SERIE_Z   ! record temporal data

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
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_LES
