!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_PARAM_LIMA_COLD
!     ###########################
!> @file
!!****  *MODD_PARAM_LIMA_COLD* - declaration of some descriptive parameters and
!!                               microphysical factors extensively used in
!!                               the LIMA cold scheme.
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13
!!      C. Barthe            14/03/2022  add CIBU and RDSF
!       J. Wurtz                03/2022: new snow characteristics
!       M. Taufour              07/2022: add concentration for snow, graupel, hail
!!
!-------------------------------------------------------------------------------
USE MODD_PARAMETERS, ONLY: JPSVNAMELGTMAX
!
IMPLICIT NONE
!
!*       1.   DESCRIPTIVE PARAMETERS
!             ----------------------
!
!     Declaration of microphysical constants, including the descriptive
!     parameters for the raindrop and the ice crystal habits, and the
!     parameters relevant of the dimensional distributions.
!
!         m(D)    = XAx * D**XBx      : Mass-MaxDim relationship
!         v(D)    = XCx * D**XDx      : Fallspeed-MaxDim relationship
!         N(Lbda) = XCCx * Lbda**XCXx : NumberConc-Slopeparam relationship
!         XF0x, XF1x, XF2x            : Ventilation factors
!         XC1x                        : Shape parameter for deposition
!
!              and
!
!         XALPHAx, XNUx                        : Generalized GAMMA law
!         Lbda = XLBx * (r_x*rho_dref)**XLBEXx : Slope parameter of the
!                                                distribution law
!
TYPE PARAM_LIMA_COLD_T
REAL      :: XLBEXI,XLBI              ! Prist. ice     distribution parameters
REAL      :: XLBEXS,XLBS,XNS          ! Snow/agg.      distribution parameters
!
REAL      :: XAI,XBI,XC_I,XDI         ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL      ::                           XF0IS,XF1IS    ! (large Di vent. coef.)
REAL      :: XDELTAI, XGAMMAI         ! cloud ice area-diameter parameters
REAL      :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
!
REAL      :: XLBDAS_MIN, XLBDAS_MAX   ! Max values allowed for the shape parameter of snow
REAL      :: XFVELOS                  ! Wurtz - snow fall speed parameterizaed after Thompson 2008
REAL      :: XTRANS_MP_GAMMAS         ! Wurtz - change between lambda value for MP and gen. gamma
!
REAL      :: XREFFI                   ! constante for ice crystal effective radius for ecRad
!
REAL, DIMENSION(:), ALLOCATABLE :: XLBEXI_SHAPE, & ! pristine ice distrib. parameters
                                   XLBI_SHAPE      ! when several shapes
REAL, DIMENSION(:), ALLOCATABLE :: XAI_SHAPE,  XBI_SHAPE, & ! a and b parameters of the mass-diam relationship
                                   XC_I_SHAPE, XDI_SHAPE, & ! c and d parameters of the fallspeed-diam relationship
                                   XGAMMAI_SHAPE, XDELTAI_SHAPE, & ! g and d parameters of proj area-diam relationship
                                   XC1I_SHAPE      ! 
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS
!             ---------------------
!
REAL      :: XFSEDRI,XFSEDCI,                  & ! Constants for sedimentation
             XFSEDRS,XFSEDCS,                  & !
             XFSEDS, XEXSEDS                     ! fluxes of ice and snow
!
REAL      :: XNUC_DEP,XEXSI_DEP,XEX_DEP,       & ! Constants for heterogeneous
             XNUC_CON,XEXTT_CON,XEX_CON,       & ! ice nucleation : DEP et CON
             XMNU0                               ! mass of nucleated ice crystal
!
REAL      :: XRHOI_HONH,XCEXP_DIFVAP_HONH,     & ! Constants for homogeneous
             XCOEF_DIFVAP_HONH,XRCOEF_HONH,    & ! haze freezing : HHONI
             XCRITSAT1_HONH,XCRITSAT2_HONH,    &
             XTMIN_HONH,XTMAX_HONH,            &
             XDLNJODT1_HONH,XDLNJODT2_HONH,    &
             XC1_HONH,XC2_HONH,XC3_HONH
!
REAL      :: XC_HONC,XR_HONC,                  & ! Constants for homogeneous
             XTEXP1_HONC,XTEXP2_HONC,          & ! droplet freezing : CHONI
             XTEXP3_HONC,XTEXP4_HONC,          &
             XTEXP5_HONC
!
REAL      :: XCSCNVI_MAX, XLBDASCNVI_MAX,      &
             XRHORSMIN,                        &
             XDSCNVI_LIM, XLBDASCNVI_LIM,      & ! Constants for snow
             XC0DEPSI,XC1DEPSI,                & ! sublimation conversion to
             XR0DEPSI,XR1DEPSI                   ! pristine ice : SCNVI
!
REAL      :: XSCFAC,                           & ! Constants for the Bergeron
             X0DEPI,X2DEPI,                    & ! Findeisen process and
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS     ! deposition
!
REAL      :: XDICNVS_LIM, XLBDAICNVS_LIM,      & ! Constants for pristine ice
             XC0DEPIS,XC1DEPIS,                & ! deposition conversion to
             XR0DEPIS,XR1DEPIS                   ! snow : ICNVS
!
REAL      :: XCOLEXIS,                         & ! Constants for snow
             XAGGS_CLARGE1,XAGGS_CLARGE2,      & ! aggregation : AGG
             XAGGS_RLARGE1,XAGGS_RLARGE2,      &
             XFIAGGS,XEXIAGGS
!
REAL      :: XACCS1, XSPONBUDS1, XSPONBUDS2,   & ! Constant for snow
             XSPONBUDS3, XSPONCOEFS2             ! spontaneous break-up
!
!??????????????????
REAL      :: XKER_ZRNIC_A1,XKER_ZRNIC_A2         ! Long-Zrnic Kernels (ini_ice_coma)
!
REAL,DIMENSION(:,:), ALLOCATABLE :: XKER_N_SSCS
REAL      :: XCOLSS,XCOLEXSS,XFNSSCS,          & !
             XLBNSSCS1,XLBNSSCS2,              & ! Constants for snow self collection
             XSCINTP1S,XSCINTP2S                 !
INTEGER      :: NSCLBDAS                         !

! for ice self-collection
REAL,DIMENSION(:,:), ALLOCATABLE :: XKER_N_ISCI
REAL      :: XSELFI, XCOLII,XCOLEXII,XFNISCI,  & !
             XLBNISCI1,XLBNISCI2,              & ! Constants for ice self collection
             XSCINTP1I,XSCINTP2I,              & !
             XSCLBDAI_MIN, XSCLBDAI_MAX
REAL, DIMENSION(:),       ALLOCATABLE :: XFISCS, XLBISCS1, XLBISCS2, XLBISCS3
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XKER_R_ISCI_SHAPE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XKER_N_ISCI_SHAPE
INTEGER   :: NSCLBDAI
!
REAL      :: XAUTO3, XAUTO4,                   & ! Constants for pristine ice
             XLAUTS,   XLAUTS_THRESHOLD,       & ! autoconversion : AUT
             XITAUTS, XITAUTS_THRESHOLD,       & ! (ini_ice_com)
             XTEXAUTI
!
REAL      :: XCONCI_MAX                          ! Limitation of the pristine
                                   ! ice concentration (init and grid-nesting)
REAL      :: XFREFFI  ! Factor to compute the cloud ice effective radius
!
! For ICE4 nucleation
REAL      :: XALPHA1
REAL      :: XALPHA2
REAL      :: XBETA1
REAL      :: XBETA2
REAL      :: XNU10
REAL      :: XNU20
! For ice crystal shapes
! pour le moment, on les met en allocatable. Peut-etre faudra t-il le modifier par la suite
! Constants for sedimentation
REAL, DIMENSION(:), ALLOCATABLE :: XFSEDCI_SHAPE      ! Constant per habit: N
REAL, DIMENSION(:), ALLOCATABLE :: XFSEDRI_SHAPE      ! Constant per habit: mr
REAL, DIMENSION(:), ALLOCATABLE :: XFSEDRI_TOT_SHAPE  ! Constant per habit: mr
! Constants for sublimation conversion to pristine ice
REAL, DIMENSION(:), ALLOCATABLE :: XC0DEPSI_SHAPE, XC1DEPSI_SHAPE, & ! sublimation conversion to
                                   XR0DEPSI_SHAPE, XR1DEPSI_SHAPE    ! pristine ice : SCNVI
! Constants for vapor deposition on ice
REAL, DIMENSION(:), ALLOCATABLE :: X0DEPI_SHAPE, & ! vapor deposition on ice
                                   X2DEPI_SHAPE
! Constants for sublimation conversion to pristine ice 
REAL, DIMENSION(:), ALLOCATABLE :: XC0DEPIS_SHAPE, XC1DEPIS_SHAPE, & ! sublimation conversion to
                                   XR0DEPIS_SHAPE, XR1DEPIS_SHAPE    ! pristine ice : SCNVI
! Constants for snow
REAL, DIMENSION(:), ALLOCATABLE :: XAGGS_RLARGE1_SHAPE, & ! Constants for snow 
                                  XAGGS_RLARGE2_SHAPE
!
REAL, DIMENSION(:), ALLOCATABLE :: XFREFFI_SHAPE  ! Factor to compute the cloud ice effective radius
!
END TYPE PARAM_LIMA_COLD_T
!
TYPE(PARAM_LIMA_COLD_T), TARGET, SAVE :: PARAM_LIMA_COLD
!
REAL, POINTER :: XLBEXI => NULL(), &
                 XLBI => NULL(), &
                 XLBEXS => NULL(), &
                 XLBS => NULL(), &
                 XNS => NULL(), &
                 XAI => NULL(), &
                 XBI => NULL(), &
                 XGAMMAI => NULL(), &
                 XDELTAI => NULL(), &
                 XC_I => NULL(), &
                 XDI => NULL(), &
                 XF0I => NULL(), &
                 XF2I => NULL(), &
                 XC1I => NULL(), &
                 XF0IS => NULL(), &
                 XF1IS => NULL(), &
                 XAS => NULL(), &
                 XBS => NULL(), &
                 XCS => NULL(), &
                 XDS => NULL(), &
                 XCCS => NULL(), &
                 XCXS => NULL(), &
                 XF0S => NULL(), &
                 XF1S => NULL(), &
                 XC1S => NULL(), &
                 XLBDAS_MIN => NULL(), &
                 XLBDAS_MAX => NULL(), &
                 XFVELOS => NULL(), &
                 XTRANS_MP_GAMMAS => NULL(), &
                 XREFFI => NULL(), &
                 XFSEDRI => NULL(), &
                 XFSEDCI => NULL(), &
                 XFSEDRS => NULL(), &
                 XFSEDCS => NULL(), &
                 XFSEDS => NULL(), &
                 XEXSEDS => NULL(), &
                 XNUC_DEP => NULL(), &
                 XEXSI_DEP => NULL(), &
                 XEX_DEP => NULL(), &
                 XNUC_CON => NULL(), &
                 XEXTT_CON => NULL(), &
                 XEX_CON => NULL(), &
                 XMNU0 => NULL(), &
                 XRHOI_HONH => NULL(), &
                 XCEXP_DIFVAP_HONH => NULL(), &
                 XCOEF_DIFVAP_HONH => NULL(), &
                 XRCOEF_HONH => NULL(), &
                 XCRITSAT1_HONH => NULL(), &
                 XCRITSAT2_HONH => NULL(), &
                 XTMIN_HONH => NULL(), &
                 XTMAX_HONH => NULL(), &
                 XDLNJODT1_HONH => NULL(), &
                 XDLNJODT2_HONH => NULL(), &
                 XC1_HONH => NULL(), &
                 XC2_HONH => NULL(), &
                 XC3_HONH => NULL(), &
                 XC_HONC => NULL(), &
                 XR_HONC => NULL(), &
                 XTEXP1_HONC => NULL(), &
                 XTEXP2_HONC => NULL(), &
                 XTEXP3_HONC => NULL(), &
                 XTEXP4_HONC => NULL(), &
                 XTEXP5_HONC => NULL(), &
                 XCSCNVI_MAX => NULL(), &
                 XLBDASCNVI_MAX => NULL(), &
                 XRHORSMIN => NULL(), &
                 XDSCNVI_LIM => NULL(), &
                 XLBDASCNVI_LIM => NULL(), &
                 XC0DEPSI => NULL(), &
                 XC1DEPSI => NULL(), &
                 XR0DEPSI => NULL(), &
                 XR1DEPSI => NULL(), &
                 XSCFAC => NULL(), &
                 X0DEPI => NULL(), &
                 X2DEPI => NULL(), &
                 X0DEPS => NULL(), &
                 X1DEPS => NULL(), &
                 XEX0DEPS => NULL(), &
                 XEX1DEPS => NULL(), &
                 XDICNVS_LIM => NULL(), &
                 XLBDAICNVS_LIM => NULL(), &
                 XC0DEPIS => NULL(), &
                 XC1DEPIS => NULL(), &
                 XR0DEPIS => NULL(), &
                 XR1DEPIS => NULL(), &
                 XCOLEXIS => NULL(), &
                 XAGGS_CLARGE1 => NULL(), &
                 XAGGS_CLARGE2 => NULL(), &
                 XAGGS_RLARGE1 => NULL(), &
                 XAGGS_RLARGE2 => NULL(), &
                 XFIAGGS => NULL(), &
                 XEXIAGGS => NULL(), &
                 XACCS1 => NULL(), &
                 XSPONBUDS1 => NULL(), &
                 XSPONBUDS2 => NULL(), &
                 XSPONBUDS3 => NULL(), &
                 XSPONCOEFS2 => NULL(), &
                 XKER_ZRNIC_A1 => NULL(), &
                 XKER_ZRNIC_A2 => NULL(), &
                 XSELFI => NULL(), &
                 XCOLSS => NULL(), &
                 XCOLEXSS => NULL(), &
                 XFNSSCS => NULL(), &
                 XLBNSSCS1 => NULL(), &
                 XLBNSSCS2 => NULL(), &
                 XSCINTP1S => NULL(), &
                 XSCINTP2S  => NULL(), &
                 XCOLII => NULL(), &
                 XCOLEXII => NULL(), &
                 XFNISCI => NULL(), &
                 XLBNISCI1 => NULL(), &
                 XLBNISCI2 => NULL(), &
                 XSCINTP1I => NULL(), &
                 XSCINTP2I  => NULL(), &
                 XSCLBDAI_MIN => NULL(),&
                 XSCLBDAI_MAX => NULL(),&                
                 XAUTO3 => NULL(), &
                 XAUTO4 => NULL(), &
                 XLAUTS => NULL(), &
                 XLAUTS_THRESHOLD => NULL(), &
                 XITAUTS => NULL(), &
                 XITAUTS_THRESHOLD => NULL(), &
                 XTEXAUTI => NULL(), &
                 XCONCI_MAX => NULL(), &
                 XFREFFI => NULL(), &
                 XALPHA1 => NULL(), &
                 XALPHA2 => NULL(), &
                 XBETA1 => NULL(), &
                 XBETA2 => NULL(), &
                 XNU10 => NULL(), &
                 XNU20 => NULL()
INTEGER, POINTER :: NSCLBDAS => NULL()
REAL,DIMENSION(:,:),POINTER :: XKER_N_SSCS => NULL()
! for ice self-collection
INTEGER, POINTER :: NSCLBDAI => NULL()
REAL,DIMENSION(:,:),    POINTER :: XKER_N_ISCI => NULL()
REAL,DIMENSION(:),      POINTER :: XLBISCS1 => NULL(), &
                                   XLBISCS2 => NULL(), &
                                   XLBISCS3 => NULL(), &
                                   XFISCS => NULL()
REAL,DIMENSION(:,:,:),  POINTER :: XKER_R_ISCI_SHAPE => NULL()
REAL,DIMENSION(:,:,:,:),POINTER :: XKER_N_ISCI_SHAPE => NULL()
! For ice crystal shapes
REAL, DIMENSION(:), POINTER :: XLBEXI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XLBI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XAI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XBI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC_I_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XDI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XGAMMAI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XDELTAI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC1I_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XFSEDCI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XFSEDRI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XFSEDRI_TOT_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC0DEPSI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC1DEPSI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XR0DEPSI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XR1DEPSI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: X0DEPI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: X2DEPI_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC0DEPIS_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XC1DEPIS_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XR0DEPIS_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XR1DEPIS_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XAGGS_RLARGE1_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XAGGS_RLARGE2_SHAPE => NULL()
REAL, DIMENSION(:), POINTER :: XFREFFI_SHAPE => NULL()
!--cb--
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(8),PARAMETER &
                              :: CLIMA_COLD_NAMES=(/'CICE    ','CSNOW   ','CGRAUPEL','CHAIL   ',&
                                                        'CIFNFREE','CIFNNUCL', &
                                                        'CCNINIMM','CCCNNUCL'/)
                                 ! basenames of the SV articles stored
                                 ! in the binary files
                                 !with IF:Ice-nuclei Free (nonactivated IFN by Dep/Cond)
                                 !     IN:Ice-nuclei Nucleated (activated IFN by Dep/Cond)
                                 !     NI:Nuclei Immersed (activated IFN by Imm)
                                 !     HF:Homogeneous Freezing
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(8),PARAMETER &
                              :: CLIMA_COLD_CONC=(/'NI ','NS ','NG ','NH ','NIF','NIN','NNI','NNH'/)!for DIAG
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(4),PARAMETER &
                              :: CLIMA_SHAPE_NAMES=(/'CICE    ','CICE2   ','CICE3   ','CICE4   '/)
                                 ! basenames of the SV articles stored
                                 ! in the binary files
                                 !with 1 - PLA : plates
                                 !     2 - COL : Columns
                                 !     3 - BUR : bullet-rosettes
                                 !     4 - POL : droxtals
CHARACTER(LEN=5),DIMENSION(4),PARAMETER &
                              :: CLIMA_SHAPE_CONC=(/'N_001','N_002','N_003','N_004'/)!for DIAG

!
CONTAINS
SUBROUTINE PARAM_LIMA_COLD_ASSOCIATE()
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('PARAM_LIMA_COLD_ASSOCIATE', 0, ZHOOK_HANDLE)
IF(.NOT. ASSOCIATED(XLBEXI)) THEN
  XLBEXI             => PARAM_LIMA_COLD%XLBEXI
  XLBI               => PARAM_LIMA_COLD%XLBI
  XLBEXS             => PARAM_LIMA_COLD%XLBEXS
  XLBS               => PARAM_LIMA_COLD%XLBS
  XNS                => PARAM_LIMA_COLD%XNS
  XAI                => PARAM_LIMA_COLD%XAI
  XBI                => PARAM_LIMA_COLD%XBI
  XGAMMAI            => PARAM_LIMA_COLD%XGAMMAI
  XDELTAI            => PARAM_LIMA_COLD%XDELTAI
  XC_I               => PARAM_LIMA_COLD%XC_I
  XDI                => PARAM_LIMA_COLD%XDI
  XF0I               => PARAM_LIMA_COLD%XF0I
  XF2I               => PARAM_LIMA_COLD%XF2I
  XC1I               => PARAM_LIMA_COLD%XC1I
  XF0IS              => PARAM_LIMA_COLD%XF0IS
  XF1IS              => PARAM_LIMA_COLD%XF1IS
  XAS                => PARAM_LIMA_COLD%XAS
  XBS                => PARAM_LIMA_COLD%XBS
  XCS                => PARAM_LIMA_COLD%XCS
  XDS                => PARAM_LIMA_COLD%XDS
  XCCS               => PARAM_LIMA_COLD%XCCS
  XCXS               => PARAM_LIMA_COLD%XCXS
  XF0S               => PARAM_LIMA_COLD%XF0S
  XF1S               => PARAM_LIMA_COLD%XF1S
  XC1S               => PARAM_LIMA_COLD%XC1S
  XLBDAS_MIN         => PARAM_LIMA_COLD%XLBDAS_MIN
  XLBDAS_MAX         => PARAM_LIMA_COLD%XLBDAS_MAX
  XFVELOS            => PARAM_LIMA_COLD%XFVELOS
  XTRANS_MP_GAMMAS   => PARAM_LIMA_COLD%XTRANS_MP_GAMMAS
  XREFFI             => PARAM_LIMA_COLD%XREFFI
  XFSEDRI            => PARAM_LIMA_COLD%XFSEDRI
  XFSEDCI            => PARAM_LIMA_COLD%XFSEDCI
  XFSEDRS            => PARAM_LIMA_COLD%XFSEDRS
  XFSEDCS            => PARAM_LIMA_COLD%XFSEDCS
  XFSEDS             => PARAM_LIMA_COLD%XFSEDS
  XEXSEDS            => PARAM_LIMA_COLD%XEXSEDS
  XNUC_DEP           => PARAM_LIMA_COLD%XNUC_DEP
  XEXSI_DEP          => PARAM_LIMA_COLD%XEXSI_DEP
  XEX_DEP            => PARAM_LIMA_COLD%XEX_DEP
  XNUC_CON           => PARAM_LIMA_COLD%XNUC_CON
  XEXTT_CON          => PARAM_LIMA_COLD%XEXTT_CON
  XEX_CON            => PARAM_LIMA_COLD%XEX_CON
  XMNU0              => PARAM_LIMA_COLD%XMNU0
  XRHOI_HONH         => PARAM_LIMA_COLD%XRHOI_HONH
  XCEXP_DIFVAP_HONH  => PARAM_LIMA_COLD%XCEXP_DIFVAP_HONH
  XCOEF_DIFVAP_HONH  => PARAM_LIMA_COLD%XCOEF_DIFVAP_HONH
  XRCOEF_HONH        => PARAM_LIMA_COLD%XRCOEF_HONH
  XCRITSAT1_HONH     => PARAM_LIMA_COLD%XCRITSAT1_HONH
  XCRITSAT2_HONH     => PARAM_LIMA_COLD%XCRITSAT2_HONH
  XTMIN_HONH         => PARAM_LIMA_COLD%XTMIN_HONH
  XTMAX_HONH         => PARAM_LIMA_COLD%XTMAX_HONH
  XDLNJODT1_HONH     => PARAM_LIMA_COLD%XDLNJODT1_HONH
  XDLNJODT2_HONH     => PARAM_LIMA_COLD%XDLNJODT2_HONH
  XC1_HONH           => PARAM_LIMA_COLD%XC1_HONH
  XC2_HONH           => PARAM_LIMA_COLD%XC2_HONH
  XC3_HONH           => PARAM_LIMA_COLD%XC3_HONH
  XC_HONC            => PARAM_LIMA_COLD%XC_HONC
  XR_HONC            => PARAM_LIMA_COLD%XR_HONC
  XTEXP1_HONC        => PARAM_LIMA_COLD%XTEXP1_HONC
  XTEXP2_HONC        => PARAM_LIMA_COLD%XTEXP2_HONC
  XTEXP3_HONC        => PARAM_LIMA_COLD%XTEXP3_HONC
  XTEXP4_HONC        => PARAM_LIMA_COLD%XTEXP4_HONC
  XTEXP5_HONC        => PARAM_LIMA_COLD%XTEXP5_HONC
  XCSCNVI_MAX        => PARAM_LIMA_COLD%XCSCNVI_MAX
  XLBDASCNVI_MAX     => PARAM_LIMA_COLD%XLBDASCNVI_MAX
  XRHORSMIN          => PARAM_LIMA_COLD%XRHORSMIN
  XDSCNVI_LIM        => PARAM_LIMA_COLD%XDSCNVI_LIM
  XLBDASCNVI_LIM     => PARAM_LIMA_COLD%XLBDASCNVI_LIM
  XC0DEPSI           => PARAM_LIMA_COLD%XC0DEPSI
  XC1DEPSI           => PARAM_LIMA_COLD%XC1DEPSI
  XR0DEPSI           => PARAM_LIMA_COLD%XR0DEPSI
  XR1DEPSI           => PARAM_LIMA_COLD%XR1DEPSI
  XSCFAC             => PARAM_LIMA_COLD%XSCFAC
  X0DEPI             => PARAM_LIMA_COLD%X0DEPI
  X2DEPI             => PARAM_LIMA_COLD%X2DEPI
  X0DEPS             => PARAM_LIMA_COLD%X0DEPS
  X1DEPS             => PARAM_LIMA_COLD%X1DEPS
  XEX0DEPS           => PARAM_LIMA_COLD%XEX0DEPS
  XEX1DEPS           => PARAM_LIMA_COLD%XEX1DEPS
  XDICNVS_LIM        => PARAM_LIMA_COLD%XDICNVS_LIM
  XLBDAICNVS_LIM     => PARAM_LIMA_COLD%XLBDAICNVS_LIM
  XC0DEPIS           => PARAM_LIMA_COLD%XC0DEPIS
  XC1DEPIS           => PARAM_LIMA_COLD%XC1DEPIS
  XR0DEPIS           => PARAM_LIMA_COLD%XR0DEPIS
  XR1DEPIS           => PARAM_LIMA_COLD%XR1DEPIS
  XCOLEXIS           => PARAM_LIMA_COLD%XCOLEXIS
  XAGGS_CLARGE1      => PARAM_LIMA_COLD%XAGGS_CLARGE1
  XAGGS_CLARGE2      => PARAM_LIMA_COLD%XAGGS_CLARGE2
  XAGGS_RLARGE1      => PARAM_LIMA_COLD%XAGGS_RLARGE1
  XAGGS_RLARGE2      => PARAM_LIMA_COLD%XAGGS_RLARGE2
  XFIAGGS            => PARAM_LIMA_COLD%XFIAGGS
  XEXIAGGS           => PARAM_LIMA_COLD%XEXIAGGS
  XACCS1             => PARAM_LIMA_COLD%XACCS1
  XSPONBUDS1         => PARAM_LIMA_COLD%XSPONBUDS1
  XSPONBUDS2         => PARAM_LIMA_COLD%XSPONBUDS2
  XSPONBUDS3         => PARAM_LIMA_COLD%XSPONBUDS3
  XSPONCOEFS2        => PARAM_LIMA_COLD%XSPONCOEFS2
  XKER_ZRNIC_A1      => PARAM_LIMA_COLD%XKER_ZRNIC_A1
  XKER_ZRNIC_A2      => PARAM_LIMA_COLD%XKER_ZRNIC_A2
  XSELFI             => PARAM_LIMA_COLD%XSELFI
  XCOLSS             => PARAM_LIMA_COLD%XCOLSS
  XCOLEXSS           => PARAM_LIMA_COLD%XCOLEXSS
  XFNSSCS            => PARAM_LIMA_COLD%XFNSSCS
  XLBNSSCS1          => PARAM_LIMA_COLD%XLBNSSCS1
  XLBNSSCS2          => PARAM_LIMA_COLD%XLBNSSCS2
  XSCINTP1S          => PARAM_LIMA_COLD%XSCINTP1S
  XSCINTP2S          => PARAM_LIMA_COLD%XSCINTP2S
  XAUTO3             => PARAM_LIMA_COLD%XAUTO3
  XAUTO4             => PARAM_LIMA_COLD%XAUTO4
  XLAUTS             => PARAM_LIMA_COLD%XLAUTS
  XLAUTS_THRESHOLD   => PARAM_LIMA_COLD%XLAUTS_THRESHOLD
  XITAUTS            => PARAM_LIMA_COLD%XITAUTS
  XITAUTS_THRESHOLD  => PARAM_LIMA_COLD%XITAUTS_THRESHOLD
  XTEXAUTI           => PARAM_LIMA_COLD%XTEXAUTI
  XCONCI_MAX         => PARAM_LIMA_COLD%XCONCI_MAX
  XFREFFI            => PARAM_LIMA_COLD%XFREFFI
  XALPHA1            => PARAM_LIMA_COLD%XALPHA1
  XALPHA2            => PARAM_LIMA_COLD%XALPHA2
  XBETA1             => PARAM_LIMA_COLD%XBETA1
  XBETA2             => PARAM_LIMA_COLD%XBETA2
  XNU10              => PARAM_LIMA_COLD%XNU10
  XNU20              => PARAM_LIMA_COLD%XNU20
  NSCLBDAS           => PARAM_LIMA_COLD%NSCLBDAS
  XCOLII             => PARAM_LIMA_COLD%XCOLII
  XCOLEXII           => PARAM_LIMA_COLD%XCOLEXII
  XFNISCI            => PARAM_LIMA_COLD%XFNISCI
  XLBNISCI1          => PARAM_LIMA_COLD%XLBNISCI1
  XLBNISCI2          => PARAM_LIMA_COLD%XLBNISCI2
  NSCLBDAI           => PARAM_LIMA_COLD%NSCLBDAI
  XSCLBDAI_MIN       => PARAM_LIMA_COLD%XSCLBDAI_MIN
  XSCLBDAI_MAX       => PARAM_LIMA_COLD%XSCLBDAI_MAX
  XSCINTP1I          => PARAM_LIMA_COLD%XSCINTP1I
  XSCINTP2I          => PARAM_LIMA_COLD%XSCINTP2I
ENDIF
IF (LHOOK) CALL DR_HOOK('PARAM_LIMA_COLD_ASSOCIATE', 1, ZHOOK_HANDLE)
END SUBROUTINE PARAM_LIMA_COLD_ASSOCIATE
!
SUBROUTINE PARAM_LIMA_COLD_ALLOCATE(HNAME, KDIM1, KDIM2)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER, INTENT(IN)          :: KDIM2

  IF (LHOOK) CALL DR_HOOK('PARAM_LIMA_COLD_ALLOCATE', 0, ZHOOK_HANDLE)
  SELECT CASE(TRIM(HNAME))
    CASE('XKER_N_SSCS')
      ALLOCATE(PARAM_LIMA_COLD%XKER_N_SSCS(KDIM1, KDIM2))
      XKER_N_SSCS => PARAM_LIMA_COLD%XKER_N_SSCS
  END SELECT
  SELECT CASE(TRIM(HNAME))
    CASE('XKER_N_ISCI')
      ALLOCATE(PARAM_LIMA_COLD%XKER_N_ISCI(KDIM1, KDIM2))
      XKER_N_ISCI => PARAM_LIMA_COLD%XKER_N_ISCI
  END SELECT
IF (LHOOK) CALL DR_HOOK('PARAM_LIMA_COLD_ALLOCATE', 1, ZHOOK_HANDLE)
END SUBROUTINE PARAM_LIMA_COLD_ALLOCATE
!
! 24/01/24
SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_1D(HNAME, KDIM1)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1

  SELECT CASE(TRIM(HNAME))
    CASE('XLBEXI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XLBEXI_SHAPE(KDIM1))
      XLBEXI_SHAPE => PARAM_LIMA_COLD%XLBEXI_SHAPE
    CASE('XLBI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XLBI_SHAPE(KDIM1))
      XLBI_SHAPE => PARAM_LIMA_COLD%XLBI_SHAPE
    CASE('XAI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XAI_SHAPE(KDIM1))
      XAI_SHAPE => PARAM_LIMA_COLD%XAI_SHAPE
    CASE('XBI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XBI_SHAPE(KDIM1))
      XBI_SHAPE => PARAM_LIMA_COLD%XBI_SHAPE
    CASE('XC_I_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC_I_SHAPE(KDIM1))
      XC_I_SHAPE => PARAM_LIMA_COLD%XC_I_SHAPE
    CASE('XDI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XDI_SHAPE(KDIM1))
      XDI_SHAPE => PARAM_LIMA_COLD%XDI_SHAPE
    CASE('XGAMMAI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XGAMMAI_SHAPE(KDIM1))
      XGAMMAI_SHAPE => PARAM_LIMA_COLD%XGAMMAI_SHAPE
    CASE('XDELTAI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XDELTAI_SHAPE(KDIM1))
      XDELTAI_SHAPE => PARAM_LIMA_COLD%XDELTAI_SHAPE      
    CASE('XC1I_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC1I_SHAPE(KDIM1))
      XC1I_SHAPE => PARAM_LIMA_COLD%XC1I_SHAPE
    CASE('XFSEDCI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XFSEDCI_SHAPE(KDIM1))
      XFSEDCI_SHAPE => PARAM_LIMA_COLD%XFSEDCI_SHAPE
    CASE('XFSEDRI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XFSEDRI_SHAPE(KDIM1))
      XFSEDRI_SHAPE => PARAM_LIMA_COLD%XFSEDRI_SHAPE
    CASE('XFSEDRI_TOT_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XFSEDRI_TOT_SHAPE(KDIM1))
      XFSEDRI_TOT_SHAPE => PARAM_LIMA_COLD%XFSEDRI_TOT_SHAPE
    CASE('XC0DEPSI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC0DEPSI_SHAPE(KDIM1))
      XC0DEPSI_SHAPE => PARAM_LIMA_COLD%XC0DEPSI_SHAPE
    CASE('XC1DEPSI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC1DEPSI_SHAPE(KDIM1))
      XC1DEPSI_SHAPE => PARAM_LIMA_COLD%XC1DEPSI_SHAPE
    CASE('XR0DEPSI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XR0DEPSI_SHAPE(KDIM1))
      XR0DEPSI_SHAPE => PARAM_LIMA_COLD%XR0DEPSI_SHAPE
    CASE('XR1DEPSI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XR1DEPSI_SHAPE(KDIM1))
      XR1DEPSI_SHAPE => PARAM_LIMA_COLD%XR1DEPSI_SHAPE
    CASE('X0DEPI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%X0DEPI_SHAPE(KDIM1))
      X0DEPI_SHAPE => PARAM_LIMA_COLD%X0DEPI_SHAPE
    CASE('X2DEPI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%X2DEPI_SHAPE(KDIM1))
      X2DEPI_SHAPE => PARAM_LIMA_COLD%X2DEPI_SHAPE
    CASE('XC0DEPIS_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC0DEPIS_SHAPE(KDIM1))
      XC0DEPIS_SHAPE => PARAM_LIMA_COLD%XC0DEPIS_SHAPE
    CASE('XC1DEPIS_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XC1DEPIS_SHAPE(KDIM1))
      XC1DEPIS_SHAPE => PARAM_LIMA_COLD%XC1DEPIS_SHAPE
    CASE('XR0DEPIS_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XR0DEPIS_SHAPE(KDIM1))
      XR0DEPIS_SHAPE => PARAM_LIMA_COLD%XR0DEPIS_SHAPE
    CASE('XR1DEPIS_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XR1DEPIS_SHAPE(KDIM1))
      XR1DEPIS_SHAPE => PARAM_LIMA_COLD%XR1DEPIS_SHAPE
    CASE('XAGGS_RLARGE1_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XAGGS_RLARGE1_SHAPE(KDIM1))
      XAGGS_RLARGE1_SHAPE => PARAM_LIMA_COLD%XAGGS_RLARGE1_SHAPE
    CASE('XAGGS_RLARGE2_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XAGGS_RLARGE2_SHAPE(KDIM1))
      XAGGS_RLARGE2_SHAPE => PARAM_LIMA_COLD%XAGGS_RLARGE2_SHAPE
    CASE('XFREFFI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XFREFFI_SHAPE(KDIM1))
      XFREFFI_SHAPE => PARAM_LIMA_COLD%XFREFFI_SHAPE
    !
    CASE('XLBISCS1')
      ALLOCATE(PARAM_LIMA_COLD%XLBISCS1(KDIM1))
      XLBISCS1 => PARAM_LIMA_COLD%XLBISCS1
    CASE('XLBISCS2')
      ALLOCATE(PARAM_LIMA_COLD%XLBISCS2(KDIM1))
      XLBISCS2 => PARAM_LIMA_COLD%XLBISCS2
    CASE('XLBISCS3')
      ALLOCATE(PARAM_LIMA_COLD%XLBISCS3(KDIM1))
      XLBISCS3 => PARAM_LIMA_COLD%XLBISCS3
    CASE('XFISCS')
      ALLOCATE(PARAM_LIMA_COLD%XFISCS(KDIM1))
      XFISCS => PARAM_LIMA_COLD%XFISCS    
    !
  END SELECT
END SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_1D
!
!
SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_3D(HNAME, KDIM1, KDIM2, KDIM3)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1
  INTEGER, INTENT(IN)          :: KDIM2
  INTEGER, INTENT(IN)          :: KDIM3  

  SELECT CASE(TRIM(HNAME))
    CASE('XKER_R_ISCI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XKER_R_ISCI_SHAPE(KDIM1, KDIM2, KDIM3))
      XKER_R_ISCI_SHAPE => PARAM_LIMA_COLD%XKER_R_ISCI_SHAPE
  END SELECT
END SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_3D
!
SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_4D(HNAME, KDIM1, KDIM2, KDIM3, KDIM4)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1
  INTEGER, INTENT(IN)          :: KDIM2
  INTEGER, INTENT(IN)          :: KDIM3  
  INTEGER, INTENT(IN)          :: KDIM4  

  SELECT CASE(TRIM(HNAME))
    CASE('XKER_N_ISCI_SHAPE')
      ALLOCATE(PARAM_LIMA_COLD%XKER_N_ISCI_SHAPE(KDIM1, KDIM2, KDIM3, KDIM4))
      XKER_N_ISCI_SHAPE => PARAM_LIMA_COLD%XKER_N_ISCI_SHAPE
  END SELECT
END SUBROUTINE PARAM_LIMA_COLD_ALLOCATE_4D
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_COLD
