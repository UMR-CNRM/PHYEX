!     ###########################
      MODULE MODD_PARAM_LIMA_COLD
!     ###########################
!
!!****  *MODD_PARAM_LIMA_COLD* - declaration of some descriptive parameters and
!!                               microphysical factors extensively used in 
!!                               the LIMA cold scheme.
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  *Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!-------------------------------------------------------------------------------
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
REAL,SAVE :: XLBEXI,XLBI              ! Prist. ice     distribution parameters
REAL,SAVE :: XLBEXS,XLBS              ! Snow/agg.      distribution parameters
!
REAL,SAVE :: XAI,XBI,XC_I,XDI         ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL,SAVE ::                           XF0IS,XF1IS    ! (large Di vent. coef.)
REAL,SAVE :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
!
REAL,SAVE :: XLBDAS_MAX               ! Max values allowed for the shape
                                      ! parameter of snow
!
CHARACTER(LEN=8),DIMENSION(5),PARAMETER &
                              :: CLIMA_COLD_NAMES=(/'CICE    ','CIFNFREE','CIFNNUCL', &
                                                        'CCNINIMM','CCCNNUCL'/)
                                 ! basenames of the SV articles stored
                                 ! in the binary files
                                 !with IF:Ice-nuclei Free (nonactivated IFN by Dep/Cond)
                                 !     IN:Ice-nuclei Nucleated (activated IFN by Dep/Cond)
                                 !     NI:Nuclei Immersed (activated IFN by Imm)
                                 !     HF:Homogeneous Freezing
CHARACTER(LEN=3),DIMENSION(5),PARAMETER &
                              :: CLIMA_COLD_CONC=(/'NI ','NIF','NIN','NNI','NNH'/)!for DIAG
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS
!             ---------------------
!
REAL,SAVE :: XFSEDRI,XFSEDCI,                  & ! Constants for sedimentation
    	     XFSEDS, XEXSEDS                     ! fluxes of ice and snow
!
REAL,SAVE :: XNUC_DEP,XEXSI_DEP,XEX_DEP,       & ! Constants for heterogeneous
    	     XNUC_CON,XEXTT_CON,XEX_CON,       & ! ice nucleation : DEP et CON
    	     XMNU0                               ! mass of nucleated ice crystal
!
REAL,SAVE :: XRHOI_HONH,XCEXP_DIFVAP_HONH,     & ! Constants for homogeneous
             XCOEF_DIFVAP_HONH,XRCOEF_HONH,    & ! haze freezing : HHONI
             XCRITSAT1_HONH,XCRITSAT2_HONH,    &
             XTMIN_HONH,XTMAX_HONH,            &
             XDLNJODT1_HONH,XDLNJODT2_HONH,    &
             XC1_HONH,XC2_HONH,XC3_HONH
!
REAL,SAVE :: XC_HONC,XR_HONC,                  & ! Constants for homogeneous
    	     XTEXP1_HONC,XTEXP2_HONC,          & ! droplet freezing : CHONI
             XTEXP3_HONC,XTEXP4_HONC,          &
             XTEXP5_HONC
!
REAL,SAVE :: XCSCNVI_MAX, XLBDASCNVI_MAX,      &
             XRHORSMIN,                        &
             XDSCNVI_LIM, XLBDASCNVI_LIM,      & ! Constants for snow
             XC0DEPSI,XC1DEPSI,                & ! sublimation conversion to
             XR0DEPSI,XR1DEPSI                   ! pristine ice : SCNVI
!
REAL,SAVE :: XSCFAC,                           & ! Constants for the Bergeron
    	     X0DEPI,X2DEPI,                    & ! Findeisen process and
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS     ! deposition
!
REAL,SAVE :: XDICNVS_LIM, XLBDAICNVS_LIM,      & ! Constants for pristine ice
             XC0DEPIS,XC1DEPIS,                & ! deposition conversion to
             XR0DEPIS,XR1DEPIS                   ! snow : ICNVS
!
REAL,SAVE :: XCOLEXIS,                         & ! Constants for snow 
    	     XAGGS_CLARGE1,XAGGS_CLARGE2,      & ! aggregation : AGG
             XAGGS_RLARGE1,XAGGS_RLARGE2
!
!??????????????????
REAL,SAVE :: XKER_ZRNIC_A1,XKER_ZRNIC_A2         ! Long-Zrnic Kernels (ini_ice_coma)
!
REAL,SAVE :: XSELFI,XCOLEXII                     ! Constants for pristine ice
                                                 ! self-collection (ini_ice_coma)
!
REAL,SAVE :: XAUTO3, XAUTO4,                   & ! Constants for pristine ice
             XLAUTS,   XLAUTS_THRESHOLD,       & ! autoconversion : AUT
             XITAUTS, XITAUTS_THRESHOLD,       & ! (ini_ice_com) 
             XTEXAUTI
!
REAL,SAVE :: XCONCI_MAX                          ! Limitation of the pristine 
                                   ! ice concentration (init and grid-nesting) 
REAL,SAVE :: XFREFFI  ! Factor to compute the cloud ice effective radius
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_COLD
