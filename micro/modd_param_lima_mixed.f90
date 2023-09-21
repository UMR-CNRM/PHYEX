!     ############################
      MODULE MODD_PARAM_LIMA_MIXED
!     ###########################{
!
!!****  *MODD_PARAM_LIMA_MIXED* - declaration of some descriptive parameters and
!!                                microphysical factors extensively used in
!!                                the LIMA mixed scheme.
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
!
IMPLICIT NONE
TYPE PARAM_LIMA_MIXED_t
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
REAL      :: XAG,XBG,XCG,XDG,XCCG,XCXG,XF0G,XF1G,XC1G ! Graupel        charact.
REAL      :: XLBEXG,XLBG,XNG          ! Graupel        distribution parameters
REAL      :: XLBDAG_MAX               ! Max values allowed for the shape
                                      ! parameter of graupeln
!
REAL      :: XAH,XBH,XCH,XDH,XCCH,XCXH,XF0H,XF1H,XC1H ! Hail           charact.
REAL      :: XALPHAH,XNUH,XLBEXH,XLBH ! Hail           distribution parameters
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS - CIBU and RDSF
!             -------------------------------------
!
! Constants for ice-ice collision : CIBU
!
REAL       :: XDCSLIM_CIBU_MIN,                & ! aggregates min diam. : 0.2 mm
              XDCSLIM_CIBU_MAX,                & ! aggregates max diam. : 1.0 mm
              XDCGLIM_CIBU_MIN,                & ! graupel min diam. : 2 mm
              XGAMINC_BOUND_CIBU_SMIN,         & ! Min val. of Lbda_s*dlim
              XGAMINC_BOUND_CIBU_SMAX,         & ! Max val. of Lbda_s*dlim
              XGAMINC_BOUND_CIBU_GMIN,         & ! Min val. of Lbda_g*dlim
              XGAMINC_BOUND_CIBU_GMAX,         & ! Max val. of Lbda_g*dlim
              XCIBUINTP_S,XCIBUINTP1_S,        & !
              XCIBUINTP2_S,                    & !
              XCIBUINTP_G,XCIBUINTP1_G,        & !
              XFACTOR_CIBU_NI,XFACTOR_CIBU_RI, & ! Factor for final CIBU Eq.
              XMOMGG_CIBU_1,XMOMGG_CIBU_2,     & ! Moment computation
              XMOMGS_CIBU_1,XMOMGS_CIBU_2,     &
              XMOMGS_CIBU_3
!
REAL, DIMENSION(:,:)      , ALLOCATABLE        &
                       :: XGAMINC_CIBU_S,      & ! Tab.incomplete Gamma function
                          XGAMINC_CIBU_G         ! Tab.incomplete Gamma function
!
! Constants for raindrop shattering : RDSF
!
REAL       :: XDCRLIM_RDSF_MIN,                & ! Raindrops min diam. : 0.2 mm
              XGAMINC_BOUND_RDSF_RMIN,         & ! Min val. of Lbda_r*dlim
              XGAMINC_BOUND_RDSF_RMAX,         & ! Max val. of Lbda_r*dlim
              XRDSFINTP_R,XRDSFINTP1_R,        & !
              XFACTOR_RDSF_NI,                 & ! Factor for final RDSF Eq.
              XMOMGR_RDSF
!
REAL, DIMENSION(:)      , ALLOCATABLE          &
                       :: XGAMINC_RDSF_R         ! Tab.incomplete Gamma function
!
!
!*       3.   MICROPHYSICAL FACTORS - Graupel
!             -------------------------------
!
REAL      :: XFSEDG, XEXSEDG, XFSEDRG, XFSEDCG   ! Sedimentation fluxes of Graupel
!
REAL      :: X0DEPG,X1DEPG,XEX0DEPG,XEX1DEPG     ! Deposition on graupel
!
REAL      :: XHMTMIN,XHMTMAX,XHM1,XHM2,        & ! Constants for the
             XHM_YIELD,XHM_COLLCS,XHM_FACTS,   & !      revised
                       XHM_COLLCG,XHM_FACTG,   & ! Hallett-Mossop process
             XGAMINC_HMC_BOUND_MIN,            & ! Min val. of Lbda_c for HMC
             XGAMINC_HMC_BOUND_MAX,            & ! Max val. of Lbda_c for HMC
             XHMSINTP1,XHMSINTP2,              & ! (this is no more used !)
             XHMLINTP1,XHMLINTP2
!
REAL      :: XDCSLIM,XCOLCS,                   & ! Constants for the riming of
             XEXCRIMSS,XCRIMSS,                & ! the aggregates : RIM
             XEXCRIMSG,XCRIMSG,                & !
             XEXSRIMCG,XSRIMCG,                & !
             XSRIMCG2, XSRIMCG3, XEXSRIMCG2,   & ! Murakami 1990
             XGAMINC_BOUND_MIN,                & ! Min val. of Lbda_s for RIM
             XGAMINC_BOUND_MAX,                & ! Max val. of Lbda_s for RIM
             XRIMINTP1,XRIMINTP2                 ! Csts for lin. interpol. of
                                                 ! the tab. incomplete Gamma law
INTEGER      :: NGAMINC                          ! Number of tab. Lbda_s
REAL, DIMENSION(:)      , ALLOCATABLE          &
                       :: XGAMINC_RIM1,        & ! Tab. incomplete Gamma funct.
                          XGAMINC_RIM2,        & ! for XDS+2 and for XBS
                          XGAMINC_RIM4,        & ! Murakami
                          XGAMINC_HMC            ! and for the HM process
!
REAL      :: XFRACCSS,                         & ! Constants for the accretion
             XFNRACCSS,                        & ! Constants for the accretion
             XLBRACCS1,XLBRACCS2,XLBRACCS3,    & ! raindrops onto the aggregates
             XLBNRACCS1,XLBNRACCS2,XLBNRACCS3, & ! raindrops onto the aggregates
             XFSACCRG,                         & ! ACC (processes RACCSS and
             XFNSACCRG,                        & ! ACC (processes RACCSS and
             XLBSACCR1,XLBSACCR2,XLBSACCR3,    & !                SACCRG)
             XLBNSACCR1,XLBNSACCR2,XLBNSACCR3, & !                SACCRG)
             XSCLBDAS_MIN,                     & ! Min val. of Lbda_s for ACC
             XSCLBDAS_MAX,                     & ! Max val. of Lbda_s for ACC
             XACCLBDAS_MIN,                    & ! Min val. of Lbda_s for ACC
             XACCLBDAS_MAX,                    & ! Max val. of Lbda_s for ACC
             XACCLBDAR_MIN,                    & ! Min val. of Lbda_r for ACC
             XACCLBDAR_MAX,                    & ! Max val. of Lbda_r for ACC
             XACCINTP1S,XACCINTP2S,            & ! Csts for bilin. interpol. of
             XACCINTP1R,XACCINTP2R               !   Lbda_s and Lbda_r in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
INTEGER      :: NACCLBDAS,                     & ! Number of Lbda_s values and
                NACCLBDAR                        !   of Lbda_r values in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
REAL,DIMENSION(:,:)      , ALLOCATABLE         &
                    :: XKER_RACCSS,       & ! Normalized kernel for RACCSS
                       XKER_RACCS,        & ! Normalized kernel for RACCS
                       XKER_SACCRG,       & ! Normalized kernel for SACCRG
                       XKER_N_RACCSS,     & ! Normalized kernel for RACCSS
                       XKER_N_RACCS,      & ! Normalized kernel for RACCS
                       XKER_N_SACCRG        ! Normalized kernel for SACCRG
REAL      :: XFSCVMG                             ! Melting-conversion factor of
                                                 ! the aggregates
!
REAL      :: XCOLIR,                           & ! Constants for rain contact
             XEXRCFRI,XRCFRI,                  & ! freezing : CFR
             XEXICFRR,XICFRR                     !
!
REAL      :: XFCDRYG,                          & ! Constants for the dry growth
             XCOLCG,                           & ! of the graupeln :
             XCOLIG,XCOLEXIG,XFIDRYG,          & !
             XCOLSG,XCOLEXSG,XFSDRYG,XFNSDRYG, & !             RCDRYG
             XLBSDRYG1,XLBSDRYG2,XLBSDRYG3,    & !             RIDRYG
             XLBNSDRYG1,XLBNSDRYG2,XLBNSDRYG3, & !             RIDRYG
             XFRDRYG,XFNRDRYG,                 & !             RSDRYG
             XLBRDRYG1,XLBRDRYG2,XLBRDRYG3,    & !             RRDRYG
             XLBNRDRYG1,XLBNRDRYG2,XLBNRDRYG3, & !             RRDRYG
             XDRYLBDAR_MIN,                    & ! Min val. of Lbda_r for DRY
             XDRYLBDAR_MAX,                    & ! Max val. of Lbda_r for DRY
             XDRYLBDAS_MIN,                    & ! Min val. of Lbda_s for DRY
             XDRYLBDAS_MAX,                    & ! Max val. of Lbda_s for DRY
             XDRYLBDAG_MIN,                    & ! Min val. of Lbda_g for DRY
             XDRYLBDAG_MAX,                    & ! Max val. of Lbda_g for DRY
             XDRYINTP1R,XDRYINTP2R,            & ! Csts for bilin. interpol. of
             XDRYINTP1S,XDRYINTP2S,            & ! Lbda_r, Lbda_s and Lbda_g in
             XDRYINTP1G,XDRYINTP2G               ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
INTEGER      :: NDRYLBDAR,                     & ! Number of Lbda_r,
                NDRYLBDAS,                     & !        of Lbda_s and
                NDRYLBDAG                        !        of Lbda_g values in
                                                 ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
REAL,DIMENSION(:,:)      , ALLOCATABLE         &
                         :: XKER_SDRYG,        & ! Normalized kernel for SDRYG
                            XKER_RDRYG,        & ! Normalized kernel for RDRYG
                            XKER_N_SDRYG,      & ! Normalized kernel for RDRYG
                            XKER_N_RDRYG         ! Normalized kernel for RDRYG
!
!-------------------------------------------------------------------------------
!
!*       4.   MICROPHYSICAL FACTORS - Hail
!             ----------------------------
!
REAL      :: XFSEDH,XEXSEDH,XFSEDRH,XFSEDCH      ! Constants for sedimentation
!
!
REAL      :: X0DEPH,X1DEPH,XEX0DEPH,XEX1DEPH     ! Constants for deposition
!
REAL      :: XFWETH,XFSWETH,XFNSWETH,          & ! Constants for the wet growth
             XLBSWETH1,XLBSWETH2,XLBSWETH3,    & ! of the hailstones : WET
             XLBNSWETH1,XLBNSWETH2,XLBNSWETH3, & ! of the hailstones : WET
             XFGWETH, XFNGWETH,                & !   processes RSWETH
             XLBGWETH1,XLBGWETH2,XLBGWETH3,    & !             RGWETH
             XLBNGWETH1,XLBNGWETH2,XLBNGWETH3, & !             RGWETH
             XWETLBDAS_MIN,                    & ! Min val. of Lbda_s for WET
             XWETLBDAS_MAX,                    & ! Max val. of Lbda_s for WET
             XWETLBDAG_MIN,                    & ! Min val. of Lbda_g for WET
             XWETLBDAG_MAX,                    & ! Max val. of Lbda_g for WET
             XWETLBDAH_MIN,                    & ! Min val. of Lbda_h for WET
             XWETLBDAH_MAX,                    & ! Max val. of Lbda_h for WET
             XWETINTP1S,XWETINTP2S,            & ! Csts for bilin. interpol. of
             XWETINTP1G,XWETINTP2G,            & ! Lbda_r, Lbda_s and Lbda_g in
             XWETINTP1H,XWETINTP2H               ! the XKER_SWETH and XKER_GWETH
                                                 !            tables
INTEGER      :: NWETLBDAS,                     & ! Number of Lbda_s,
                NWETLBDAG,                     & !        of Lbda_g and
                NWETLBDAH                        !        of Lbda_h values in
                                                 ! the XKER_SWETH and XKER_GWETH
                                                 !            tables
REAL,DIMENSION(:,:), ALLOCATABLE               &
                         :: XKER_SWETH,        & ! Normalized kernel for SWETH
                            XKER_GWETH,        & ! Normalized kernel for GWETH
                            XKER_N_SWETH,      & ! Normalized kernel for GWETH
                            XKER_N_GWETH         ! Normalized kernel for GWETH
END TYPE PARAM_LIMA_MIXED_t
!
TYPE(PARAM_LIMA_MIXED_t), TARGET       :: PARAM_LIMA_MIXED
!
REAL, POINTER :: XAG => NULL(), &
                 XBG => NULL(), &
                 XCG => NULL(), &
                 XDG => NULL(), &
                 XCCG => NULL(), &
                 XCXG => NULL(), &
                 XF0G => NULL(), &
                 XF1G => NULL(), &
                 XC1G => NULL(), &
                 XLBEXG => NULL(), &
                 XLBG => NULL(), &
                 XNG => NULL(), &
                 XLBDAG_MAX => NULL(), &
                 XAH => NULL(), &
                 XBH => NULL(), &
                 XCH => NULL(), &
                 XDH => NULL(), &
                 XCCH => NULL(), &
                 XCXH => NULL(), &
                 XF0H => NULL(), &
                 XF1H => NULL(), &
                 XC1H => NULL(), &
                 XALPHAH => NULL(), &
                 XNUH => NULL(), &
                 XLBEXH => NULL(), &
                 XLBH => NULL(), &
                 XDCSLIM_CIBU_MIN => NULL(), &
                 XDCSLIM_CIBU_MAX => NULL(), &
                 XDCGLIM_CIBU_MIN => NULL(), &
                 XGAMINC_BOUND_CIBU_SMIN => NULL(), &
                 XGAMINC_BOUND_CIBU_SMAX => NULL(), &
                 XGAMINC_BOUND_CIBU_GMIN => NULL(), &
                 XGAMINC_BOUND_CIBU_GMAX => NULL(), &
                 XCIBUINTP_S => NULL(), &
                 XCIBUINTP1_S => NULL(), &
                 XCIBUINTP2_S => NULL(), &
                 XCIBUINTP_G => NULL(), &
                 XCIBUINTP1_G => NULL(), &
                 XFACTOR_CIBU_NI => NULL(), &
                 XFACTOR_CIBU_RI => NULL(), &
                 XMOMGG_CIBU_1 => NULL(), &
                 XMOMGG_CIBU_2 => NULL(), &
                 XMOMGS_CIBU_1 => NULL(), &
                 XMOMGS_CIBU_2 => NULL(), &
                 XMOMGS_CIBU_3 => NULL(), &
                 XDCRLIM_RDSF_MIN => NULL(), &
                 XGAMINC_BOUND_RDSF_RMIN => NULL(), &
                 XGAMINC_BOUND_RDSF_RMAX => NULL(), &
                 XRDSFINTP_R => NULL(), &
                 XRDSFINTP1_R => NULL(), &
                 XFACTOR_RDSF_NI => NULL(), &
                 XMOMGR_RDSF => NULL(), &
                 XFSEDG => NULL(), &
                 XEXSEDG => NULL(), &
                 XFSEDRG => NULL(), &
                 XFSEDCG => NULL(), &
                 X0DEPG => NULL(), &
                 X1DEPG => NULL(), &
                 XEX0DEPG => NULL(), &
                 XEX1DEPG => NULL(), &
                 XHMTMIN => NULL(), &
                 XHMTMAX => NULL(), &
                 XHM1 => NULL(), &
                 XHM2 => NULL(), &
                 XHM_YIELD => NULL(), &
                 XHM_COLLCS => NULL(), &
                 XHM_FACTS => NULL(), &
                 XHM_COLLCG => NULL(), &
                 XHM_FACTG => NULL(), &
                 XGAMINC_HMC_BOUND_MIN => NULL(), &
                 XGAMINC_HMC_BOUND_MAX => NULL(), &
                 XHMSINTP1 => NULL(), &
                 XHMSINTP2 => NULL(), &
                 XHMLINTP1 => NULL(), &
                 XHMLINTP2 => NULL(), &
                 XDCSLIM => NULL(), &
                 XCOLCS => NULL(), &
                 XEXCRIMSS => NULL(), &
                 XCRIMSS => NULL(), &
                 XEXCRIMSG => NULL(), &
                 XCRIMSG => NULL(), &
                 XEXSRIMCG => NULL(), &
                 XSRIMCG => NULL(), &
                 XSRIMCG2 => NULL(), &
                 XSRIMCG3 => NULL(), &
                 XEXSRIMCG2 => NULL(), &
                 XGAMINC_BOUND_MIN => NULL(), &
                 XGAMINC_BOUND_MAX => NULL(), &
                 XRIMINTP1 => NULL(), &
                 XRIMINTP2 => NULL(), &
                 XFRACCSS => NULL(), &
                 XFNRACCSS => NULL(), &
                 XLBRACCS1 => NULL(), &
                 XLBRACCS2 => NULL(), &
                 XLBRACCS3 => NULL(), &
                 XLBNRACCS1 => NULL(), &
                 XLBNRACCS2 => NULL(), &
                 XLBNRACCS3 => NULL(), &
                 XFSACCRG => NULL(), &
                 XFNSACCRG => NULL(), &
                 XLBSACCR1 => NULL(), &
                 XLBSACCR2 => NULL(), &
                 XLBSACCR3 => NULL(), &
                 XLBNSACCR1 => NULL(), &
                 XLBNSACCR2 => NULL(), &
                 XLBNSACCR3 => NULL(), &
                 XSCLBDAS_MIN => NULL(), &
                 XSCLBDAS_MAX => NULL(), &
                 XACCLBDAS_MIN => NULL(), &
                 XACCLBDAS_MAX => NULL(), &
                 XACCLBDAR_MIN => NULL(), &
                 XACCLBDAR_MAX => NULL(), &
                 XACCINTP1S => NULL(), &
                 XACCINTP2S => NULL(), &
                 XACCINTP1R => NULL(), &
                 XACCINTP2R => NULL(), &
                 XFSCVMG => NULL(), &
                 XCOLIR => NULL(), &
                 XEXRCFRI => NULL(), &
                 XRCFRI => NULL(), &
                 XEXICFRR => NULL(), &
                 XICFRR => NULL(), &
                 XFCDRYG => NULL(), &
                 XCOLCG => NULL(), &
                 XCOLIG => NULL(), &
                 XCOLEXIG => NULL(), &
                 XFIDRYG => NULL(), &
                 XCOLSG => NULL(), &
                 XCOLEXSG => NULL(), &
                 XFSDRYG => NULL(), &
                 XFNSDRYG => NULL(), &
                 XLBSDRYG1 => NULL(), &
                 XLBSDRYG2 => NULL(), &
                 XLBSDRYG3 => NULL(), &
                 XLBNSDRYG1 => NULL(), &
                 XLBNSDRYG2 => NULL(), &
                 XLBNSDRYG3 => NULL(), &
                 XFRDRYG => NULL(), &
                 XFNRDRYG => NULL(), &
                 XLBRDRYG1 => NULL(), &
                 XLBRDRYG2 => NULL(), &
                 XLBRDRYG3 => NULL(), &
                 XLBNRDRYG1 => NULL(), &
                 XLBNRDRYG2 => NULL(), &
                 XLBNRDRYG3 => NULL(), &
                 XDRYLBDAR_MIN => NULL(), &
                 XDRYLBDAR_MAX => NULL(), &
                 XDRYLBDAS_MIN => NULL(), &
                 XDRYLBDAS_MAX => NULL(), &
                 XDRYLBDAG_MIN => NULL(), &
                 XDRYLBDAG_MAX => NULL(), &
                 XDRYINTP1R => NULL(), &
                 XDRYINTP2R => NULL(), &
                 XDRYINTP1S => NULL(), &
                 XDRYINTP2S => NULL(), &
                 XDRYINTP1G => NULL(), &
                 XDRYINTP2G => NULL(), &
                 XFSEDH => NULL(), &
                 XEXSEDH => NULL(), &
                 XFSEDRH => NULL(), &
                 XFSEDCH => NULL(), &
                 X0DEPH => NULL(), &
                 X1DEPH => NULL(), &
                 XEX0DEPH => NULL(), &
                 XEX1DEPH => NULL(), &
                 XFWETH => NULL(), &
                 XFSWETH => NULL(), &
                 XFNSWETH => NULL(), &
                 XLBSWETH1 => NULL(), &
                 XLBSWETH2 => NULL(), &
                 XLBSWETH3 => NULL(), &
                 XLBNSWETH1 => NULL(), &
                 XLBNSWETH2 => NULL(), &
                 XLBNSWETH3 => NULL(), &
                 XFGWETH => NULL(), &
                 XFNGWETH => NULL(), &
                 XLBGWETH1 => NULL(), &
                 XLBGWETH2 => NULL(), &
                 XLBGWETH3 => NULL(), &
                 XLBNGWETH1 => NULL(), &
                 XLBNGWETH2 => NULL(), &
                 XLBNGWETH3 => NULL(), &
                 XWETLBDAS_MIN => NULL(), &
                 XWETLBDAS_MAX => NULL(), &
                 XWETLBDAG_MIN => NULL(), &
                 XWETLBDAG_MAX => NULL(), &
                 XWETLBDAH_MIN => NULL(), &
                 XWETLBDAH_MAX => NULL(), &
                 XWETINTP1S => NULL(), &
                 XWETINTP2S => NULL(), &
                 XWETINTP1G => NULL(), &
                 XWETINTP2G => NULL(), &
                 XWETINTP1H => NULL(), &
                 XWETINTP2H => NULL()

INTEGER, POINTER :: NGAMINC => NULL(), &
                    NACCLBDAS => NULL(), &
                    NACCLBDAR => NULL(), &
                    NDRYLBDAR => NULL(), &
                    NDRYLBDAS => NULL(), &
                    NDRYLBDAG => NULL(), &
                    NWETLBDAS => NULL(), &
                    NWETLBDAG => NULL(), &
                    NWETLBDAH => NULL()

REAL, DIMENSION(:), POINTER :: XGAMINC_RDSF_R => NULL(), &
                               XGAMINC_RIM1 => NULL(), &
                               XGAMINC_RIM2 => NULL(), &
                               XGAMINC_RIM4 => NULL(), &
                               XGAMINC_HMC => NULL()
REAL, DIMENSION(:,:), POINTER :: XGAMINC_CIBU_S => NULL(), &
                                 XGAMINC_CIBU_G => NULL(), &
                                 XKER_RACCSS => NULL(), &
                                 XKER_RACCS => NULL(), &
                                 XKER_SACCRG => NULL(), &
                                 XKER_N_RACCSS => NULL(), &
                                 XKER_N_RACCS => NULL(), &
                                 XKER_N_SACCRG => NULL(), &
                                 XKER_SDRYG => NULL(), &
                                 XKER_RDRYG => NULL(), &
                                 XKER_N_SDRYG => NULL(), &
                                 XKER_N_RDRYG => NULL(), &
                                 XKER_SWETH => NULL(), &
                                 XKER_GWETH => NULL(), &
                                 XKER_N_SWETH => NULL(), &
                                 XKER_N_GWETH => NULL()
CONTAINS
SUBROUTINE PARAM_LIMA_MIXED_ASSOCIATE()
IF(.NOT. ASSOCIATED(XAG)) THEN
  XAG                      => PARAM_LIMA_MIXED%XAG
  XBG                      => PARAM_LIMA_MIXED%XBG
  XCG                      => PARAM_LIMA_MIXED%XCG
  XDG                      => PARAM_LIMA_MIXED%XDG
  XCCG                     => PARAM_LIMA_MIXED%XCCG
  XCXG                     => PARAM_LIMA_MIXED%XCXG
  XF0G                     => PARAM_LIMA_MIXED%XF0G
  XF1G                     => PARAM_LIMA_MIXED%XF1G
  XC1G                     => PARAM_LIMA_MIXED%XC1G
  XLBEXG                   => PARAM_LIMA_MIXED%XLBEXG
  XLBG                     => PARAM_LIMA_MIXED%XLBG
  XNG                      => PARAM_LIMA_MIXED%XNG
  XLBDAG_MAX               => PARAM_LIMA_MIXED%XLBDAG_MAX
  XAH                      => PARAM_LIMA_MIXED%XAH
  XBH                      => PARAM_LIMA_MIXED%XBH
  XCH                      => PARAM_LIMA_MIXED%XCH
  XDH                      => PARAM_LIMA_MIXED%XDH
  XCCH                     => PARAM_LIMA_MIXED%XCCH
  XCXH                     => PARAM_LIMA_MIXED%XCXH
  XF0H                     => PARAM_LIMA_MIXED%XF0H
  XF1H                     => PARAM_LIMA_MIXED%XF1H
  XC1H                     => PARAM_LIMA_MIXED%XC1H
  XALPHAH                  => PARAM_LIMA_MIXED%XALPHAH
  XNUH                     => PARAM_LIMA_MIXED%XNUH
  XLBEXH                   => PARAM_LIMA_MIXED%XLBEXH
  XLBH                     => PARAM_LIMA_MIXED%XLBH
  XDCSLIM_CIBU_MIN         => PARAM_LIMA_MIXED%XDCSLIM_CIBU_MIN
  XDCSLIM_CIBU_MAX         => PARAM_LIMA_MIXED%XDCSLIM_CIBU_MAX
  XDCGLIM_CIBU_MIN         => PARAM_LIMA_MIXED%XDCGLIM_CIBU_MIN
  XGAMINC_BOUND_CIBU_SMIN  => PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMIN
  XGAMINC_BOUND_CIBU_SMAX  => PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMAX
  XGAMINC_BOUND_CIBU_GMIN  => PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMIN
  XGAMINC_BOUND_CIBU_GMAX  => PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMAX
  XCIBUINTP_S              => PARAM_LIMA_MIXED%XCIBUINTP_S
  XCIBUINTP1_S             => PARAM_LIMA_MIXED%XCIBUINTP1_S
  XCIBUINTP2_S             => PARAM_LIMA_MIXED%XCIBUINTP2_S
  XCIBUINTP_G              => PARAM_LIMA_MIXED%XCIBUINTP_G
  XCIBUINTP1_G             => PARAM_LIMA_MIXED%XCIBUINTP1_G
  XFACTOR_CIBU_NI          => PARAM_LIMA_MIXED%XFACTOR_CIBU_NI
  XFACTOR_CIBU_RI          => PARAM_LIMA_MIXED%XFACTOR_CIBU_RI
  XMOMGG_CIBU_1            => PARAM_LIMA_MIXED%XMOMGG_CIBU_1
  XMOMGG_CIBU_2            => PARAM_LIMA_MIXED%XMOMGG_CIBU_2
  XMOMGS_CIBU_1            => PARAM_LIMA_MIXED%XMOMGS_CIBU_1
  XMOMGS_CIBU_2            => PARAM_LIMA_MIXED%XMOMGS_CIBU_2
  XMOMGS_CIBU_3            => PARAM_LIMA_MIXED%XMOMGS_CIBU_3
  XDCRLIM_RDSF_MIN         => PARAM_LIMA_MIXED%XDCRLIM_RDSF_MIN
  XGAMINC_BOUND_RDSF_RMIN  => PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMIN
  XGAMINC_BOUND_RDSF_RMAX  => PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMAX
  XRDSFINTP_R              => PARAM_LIMA_MIXED%XRDSFINTP_R
  XRDSFINTP1_R             => PARAM_LIMA_MIXED%XRDSFINTP1_R
  XFACTOR_RDSF_NI          => PARAM_LIMA_MIXED%XFACTOR_RDSF_NI
  XMOMGR_RDSF              => PARAM_LIMA_MIXED%XMOMGR_RDSF
  XFSEDG                   => PARAM_LIMA_MIXED%XFSEDG
  XEXSEDG                  => PARAM_LIMA_MIXED%XEXSEDG
  XFSEDRG                  => PARAM_LIMA_MIXED%XFSEDRG
  XFSEDCG                  => PARAM_LIMA_MIXED%XFSEDCG
  X0DEPG                   => PARAM_LIMA_MIXED%X0DEPG
  X1DEPG                   => PARAM_LIMA_MIXED%X1DEPG
  XEX0DEPG                 => PARAM_LIMA_MIXED%XEX0DEPG
  XEX1DEPG                 => PARAM_LIMA_MIXED%XEX1DEPG
  XHMTMIN                  => PARAM_LIMA_MIXED%XHMTMIN
  XHMTMAX                  => PARAM_LIMA_MIXED%XHMTMAX
  XHM1                     => PARAM_LIMA_MIXED%XHM1
  XHM2                     => PARAM_LIMA_MIXED%XHM2
  XHM_YIELD                => PARAM_LIMA_MIXED%XHM_YIELD
  XHM_COLLCS               => PARAM_LIMA_MIXED%XHM_COLLCS
  XHM_FACTS                => PARAM_LIMA_MIXED%XHM_FACTS
  XHM_COLLCG               => PARAM_LIMA_MIXED%XHM_COLLCG
  XHM_FACTG                => PARAM_LIMA_MIXED%XHM_FACTG
  XGAMINC_HMC_BOUND_MIN    => PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MIN
  XGAMINC_HMC_BOUND_MAX    => PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MAX
  XHMSINTP1                => PARAM_LIMA_MIXED%XHMSINTP1
  XHMSINTP2                => PARAM_LIMA_MIXED%XHMSINTP2
  XHMLINTP1                => PARAM_LIMA_MIXED%XHMLINTP1
  XHMLINTP2                => PARAM_LIMA_MIXED%XHMLINTP2
  XDCSLIM                  => PARAM_LIMA_MIXED%XDCSLIM
  XCOLCS                   => PARAM_LIMA_MIXED%XCOLCS
  XEXCRIMSS                => PARAM_LIMA_MIXED%XEXCRIMSS
  XCRIMSS                  => PARAM_LIMA_MIXED%XCRIMSS
  XEXCRIMSG                => PARAM_LIMA_MIXED%XEXCRIMSG
  XCRIMSG                  => PARAM_LIMA_MIXED%XCRIMSG
  XEXSRIMCG                => PARAM_LIMA_MIXED%XEXSRIMCG
  XSRIMCG                  => PARAM_LIMA_MIXED%XSRIMCG
  XSRIMCG2                 => PARAM_LIMA_MIXED%XSRIMCG2
  XSRIMCG3                 => PARAM_LIMA_MIXED%XSRIMCG3
  XEXSRIMCG2               => PARAM_LIMA_MIXED%XEXSRIMCG2
  XGAMINC_BOUND_MIN        => PARAM_LIMA_MIXED%XGAMINC_BOUND_MIN
  XGAMINC_BOUND_MAX        => PARAM_LIMA_MIXED%XGAMINC_BOUND_MAX
  XRIMINTP1                => PARAM_LIMA_MIXED%XRIMINTP1
  XRIMINTP2                => PARAM_LIMA_MIXED%XRIMINTP2
  XFRACCSS                 => PARAM_LIMA_MIXED%XFRACCSS
  XFNRACCSS                => PARAM_LIMA_MIXED%XFNRACCSS
  XLBRACCS1                => PARAM_LIMA_MIXED%XLBRACCS1
  XLBRACCS2                => PARAM_LIMA_MIXED%XLBRACCS2
  XLBRACCS3                => PARAM_LIMA_MIXED%XLBRACCS3
  XLBNRACCS1               => PARAM_LIMA_MIXED%XLBNRACCS1
  XLBNRACCS2               => PARAM_LIMA_MIXED%XLBNRACCS2
  XLBNRACCS3               => PARAM_LIMA_MIXED%XLBNRACCS3
  XFSACCRG                 => PARAM_LIMA_MIXED%XFSACCRG
  XFNSACCRG                => PARAM_LIMA_MIXED%XFNSACCRG
  XLBSACCR1                => PARAM_LIMA_MIXED%XLBSACCR1
  XLBSACCR2                => PARAM_LIMA_MIXED%XLBSACCR2
  XLBSACCR3                => PARAM_LIMA_MIXED%XLBSACCR3
  XLBNSACCR1               => PARAM_LIMA_MIXED%XLBNSACCR1
  XLBNSACCR2               => PARAM_LIMA_MIXED%XLBNSACCR2
  XLBNSACCR3               => PARAM_LIMA_MIXED%XLBNSACCR3
  XSCLBDAS_MIN             => PARAM_LIMA_MIXED%XSCLBDAS_MIN
  XSCLBDAS_MAX             => PARAM_LIMA_MIXED%XSCLBDAS_MAX
  XACCLBDAS_MIN            => PARAM_LIMA_MIXED%XACCLBDAS_MIN
  XACCLBDAS_MAX            => PARAM_LIMA_MIXED%XACCLBDAS_MAX
  XACCLBDAR_MIN            => PARAM_LIMA_MIXED%XACCLBDAR_MIN
  XACCLBDAR_MAX            => PARAM_LIMA_MIXED%XACCLBDAR_MAX
  XACCINTP1S               => PARAM_LIMA_MIXED%XACCINTP1S
  XACCINTP2S               => PARAM_LIMA_MIXED%XACCINTP2S
  XACCINTP1R               => PARAM_LIMA_MIXED%XACCINTP1R
  XACCINTP2R               => PARAM_LIMA_MIXED%XACCINTP2R
  XFSCVMG                  => PARAM_LIMA_MIXED%XFSCVMG
  XCOLIR                   => PARAM_LIMA_MIXED%XCOLIR
  XEXRCFRI                 => PARAM_LIMA_MIXED%XEXRCFRI
  XRCFRI                   => PARAM_LIMA_MIXED%XRCFRI
  XEXICFRR                 => PARAM_LIMA_MIXED%XEXICFRR
  XICFRR                   => PARAM_LIMA_MIXED%XICFRR
  XFCDRYG                  => PARAM_LIMA_MIXED%XFCDRYG
  XCOLCG                   => PARAM_LIMA_MIXED%XCOLCG
  XCOLIG                   => PARAM_LIMA_MIXED%XCOLIG
  XCOLEXIG                 => PARAM_LIMA_MIXED%XCOLEXIG
  XFIDRYG                  => PARAM_LIMA_MIXED%XFIDRYG
  XCOLSG                   => PARAM_LIMA_MIXED%XCOLSG
  XCOLEXSG                 => PARAM_LIMA_MIXED%XCOLEXSG
  XFSDRYG                  => PARAM_LIMA_MIXED%XFSDRYG
  XFNSDRYG                 => PARAM_LIMA_MIXED%XFNSDRYG
  XLBSDRYG1                => PARAM_LIMA_MIXED%XLBSDRYG1
  XLBSDRYG2                => PARAM_LIMA_MIXED%XLBSDRYG2
  XLBSDRYG3                => PARAM_LIMA_MIXED%XLBSDRYG3
  XLBNSDRYG1               => PARAM_LIMA_MIXED%XLBNSDRYG1
  XLBNSDRYG2               => PARAM_LIMA_MIXED%XLBNSDRYG2
  XLBNSDRYG3               => PARAM_LIMA_MIXED%XLBNSDRYG3
  XFRDRYG                  => PARAM_LIMA_MIXED%XFRDRYG
  XFNRDRYG                 => PARAM_LIMA_MIXED%XFNRDRYG
  XLBRDRYG1                => PARAM_LIMA_MIXED%XLBRDRYG1
  XLBRDRYG2                => PARAM_LIMA_MIXED%XLBRDRYG2
  XLBRDRYG3                => PARAM_LIMA_MIXED%XLBRDRYG3
  XLBNRDRYG1               => PARAM_LIMA_MIXED%XLBNRDRYG1
  XLBNRDRYG2               => PARAM_LIMA_MIXED%XLBNRDRYG2
  XLBNRDRYG3               => PARAM_LIMA_MIXED%XLBNRDRYG3
  XDRYLBDAR_MIN            => PARAM_LIMA_MIXED%XDRYLBDAR_MIN
  XDRYLBDAR_MAX            => PARAM_LIMA_MIXED%XDRYLBDAR_MAX
  XDRYLBDAS_MIN            => PARAM_LIMA_MIXED%XDRYLBDAS_MIN
  XDRYLBDAS_MAX            => PARAM_LIMA_MIXED%XDRYLBDAS_MAX
  XDRYLBDAG_MIN            => PARAM_LIMA_MIXED%XDRYLBDAG_MIN
  XDRYLBDAG_MAX            => PARAM_LIMA_MIXED%XDRYLBDAG_MAX
  XDRYINTP1R               => PARAM_LIMA_MIXED%XDRYINTP1R
  XDRYINTP2R               => PARAM_LIMA_MIXED%XDRYINTP2R
  XDRYINTP1S               => PARAM_LIMA_MIXED%XDRYINTP1S
  XDRYINTP2S               => PARAM_LIMA_MIXED%XDRYINTP2S
  XDRYINTP1G               => PARAM_LIMA_MIXED%XDRYINTP1G
  XDRYINTP2G               => PARAM_LIMA_MIXED%XDRYINTP2G
  XFSEDH                   => PARAM_LIMA_MIXED%XFSEDH
  XEXSEDH                  => PARAM_LIMA_MIXED%XEXSEDH
  XFSEDRH                  => PARAM_LIMA_MIXED%XFSEDRH
  XFSEDCH                  => PARAM_LIMA_MIXED%XFSEDCH
  X0DEPH                   => PARAM_LIMA_MIXED%X0DEPH
  X1DEPH                   => PARAM_LIMA_MIXED%X1DEPH
  XEX0DEPH                 => PARAM_LIMA_MIXED%XEX0DEPH
  XEX1DEPH                 => PARAM_LIMA_MIXED%XEX1DEPH
  XFWETH                   => PARAM_LIMA_MIXED%XFWETH
  XFSWETH                  => PARAM_LIMA_MIXED%XFSWETH
  XFNSWETH                 => PARAM_LIMA_MIXED%XFNSWETH
  XLBSWETH1                => PARAM_LIMA_MIXED%XLBSWETH1
  XLBSWETH2                => PARAM_LIMA_MIXED%XLBSWETH2
  XLBSWETH3                => PARAM_LIMA_MIXED%XLBSWETH3
  XLBNSWETH1               => PARAM_LIMA_MIXED%XLBNSWETH1
  XLBNSWETH2               => PARAM_LIMA_MIXED%XLBNSWETH2
  XLBNSWETH3               => PARAM_LIMA_MIXED%XLBNSWETH3
  XFGWETH                  => PARAM_LIMA_MIXED%XFGWETH
  XFNGWETH                 => PARAM_LIMA_MIXED%XFNGWETH
  XLBGWETH1                => PARAM_LIMA_MIXED%XLBGWETH1
  XLBGWETH2                => PARAM_LIMA_MIXED%XLBGWETH2
  XLBGWETH3                => PARAM_LIMA_MIXED%XLBGWETH3
  XLBNGWETH1               => PARAM_LIMA_MIXED%XLBNGWETH1
  XLBNGWETH2               => PARAM_LIMA_MIXED%XLBNGWETH2
  XLBNGWETH3               => PARAM_LIMA_MIXED%XLBNGWETH3
  XWETLBDAS_MIN            => PARAM_LIMA_MIXED%XWETLBDAS_MIN
  XWETLBDAS_MAX            => PARAM_LIMA_MIXED%XWETLBDAS_MAX
  XWETLBDAG_MIN            => PARAM_LIMA_MIXED%XWETLBDAG_MIN
  XWETLBDAG_MAX            => PARAM_LIMA_MIXED%XWETLBDAG_MAX
  XWETLBDAH_MIN            => PARAM_LIMA_MIXED%XWETLBDAH_MIN
  XWETLBDAH_MAX            => PARAM_LIMA_MIXED%XWETLBDAH_MAX
  XWETINTP1S               => PARAM_LIMA_MIXED%XWETINTP1S
  XWETINTP2S               => PARAM_LIMA_MIXED%XWETINTP2S
  XWETINTP1G               => PARAM_LIMA_MIXED%XWETINTP1G
  XWETINTP2G               => PARAM_LIMA_MIXED%XWETINTP2G
  XWETINTP1H               => PARAM_LIMA_MIXED%XWETINTP1H
  XWETINTP2H               => PARAM_LIMA_MIXED%XWETINTP2H

  NGAMINC                  => PARAM_LIMA_MIXED%NGAMINC
  NACCLBDAS                => PARAM_LIMA_MIXED%NACCLBDAS
  NACCLBDAR                => PARAM_LIMA_MIXED%NACCLBDAR
  NDRYLBDAR                => PARAM_LIMA_MIXED%NDRYLBDAR
  NDRYLBDAS                => PARAM_LIMA_MIXED%NDRYLBDAS
  NDRYLBDAG                => PARAM_LIMA_MIXED%NDRYLBDAG
  NWETLBDAS                => PARAM_LIMA_MIXED%NWETLBDAS
  NWETLBDAG                => PARAM_LIMA_MIXED%NWETLBDAG
  NWETLBDAH                => PARAM_LIMA_MIXED%NWETLBDAH
ENDIF
END SUBROUTINE PARAM_LIMA_MIXED_ASSOCIATE
!
SUBROUTINE PARAM_LIMA_MIXED_ALLOCATE(HNAME, KDIM1, KDIM2)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1
  INTEGER, OPTIONAL, INTENT(IN):: KDIM2

  SELECT CASE(TRIM(HNAME))
    !1d
    CASE('XGAMINC_RDSF_R')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_RDSF_R(KDIM1))
      XGAMINC_RDSF_R => PARAM_LIMA_MIXED%XGAMINC_RDSF_R
    CASE('XGAMINC_RIM1')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_RIM1(KDIM1))
      XGAMINC_RIM1 => PARAM_LIMA_MIXED%XGAMINC_RIM1
    CASE('XGAMINC_RIM2')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_RIM2(KDIM1))
      XGAMINC_RIM2 => PARAM_LIMA_MIXED%XGAMINC_RIM2
    CASE('XGAMINC_RIM4')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_RIM4(KDIM1))
      XGAMINC_RIM4 => PARAM_LIMA_MIXED%XGAMINC_RIM4
    CASE('XGAMINC_HMC')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_HMC(KDIM1))
      XGAMINC_HMC => PARAM_LIMA_MIXED%XGAMINC_HMC
   !2d
   CASE('XGAMINC_CIBU_S')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_CIBU_S(KDIM1, KDIM2))
      XGAMINC_CIBU_S => PARAM_LIMA_MIXED%XGAMINC_CIBU_S
   CASE('XGAMINC_CIBU_G')
      ALLOCATE(PARAM_LIMA_MIXED%XGAMINC_CIBU_G(KDIM1, KDIM2))
      XGAMINC_CIBU_G => PARAM_LIMA_MIXED%XGAMINC_CIBU_G
   CASE('XKER_RACCSS')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_RACCSS(KDIM1, KDIM2))
      XKER_RACCSS => PARAM_LIMA_MIXED%XKER_RACCSS
  CASE('XKER_RACCS')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_RACCS(KDIM1, KDIM2))
      XKER_RACCS => PARAM_LIMA_MIXED%XKER_RACCS
  CASE('XKER_SACCRG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_SACCRG(KDIM1, KDIM2))
      XKER_SACCRG => PARAM_LIMA_MIXED%XKER_SACCRG
  CASE('XKER_N_RACCSS')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_RACCSS(KDIM1, KDIM2))
      XKER_N_RACCSS => PARAM_LIMA_MIXED%XKER_N_RACCSS
  CASE('XKER_N_RACCS')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_RACCS(KDIM1, KDIM2))
      XKER_N_RACCS => PARAM_LIMA_MIXED%XKER_N_RACCS
  CASE('XKER_N_SACCRG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_SACCRG(KDIM1, KDIM2))
      XKER_N_SACCRG => PARAM_LIMA_MIXED%XKER_N_SACCRG
  CASE('XKER_SDRYG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_SDRYG(KDIM1, KDIM2))
      XKER_SDRYG => PARAM_LIMA_MIXED%XKER_SDRYG
  CASE('XKER_RDRYG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_RDRYG(KDIM1, KDIM2))
      XKER_RDRYG => PARAM_LIMA_MIXED%XKER_RDRYG
  CASE('XKER_N_SDRYG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_SDRYG(KDIM1, KDIM2))
      XKER_N_SDRYG => PARAM_LIMA_MIXED%XKER_N_SDRYG
  CASE('XKER_N_RDRYG')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_RDRYG(KDIM1, KDIM2))
      XKER_N_RDRYG => PARAM_LIMA_MIXED%XKER_N_RDRYG
  CASE('XKER_SWETH')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_SWETH(KDIM1, KDIM2))
      XKER_SWETH => PARAM_LIMA_MIXED%XKER_SWETH
  CASE('XKER_GWETH')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_GWETH(KDIM1, KDIM2))
      XKER_GWETH => PARAM_LIMA_MIXED%XKER_GWETH
  CASE('XKER_N_SWETH')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_SWETH(KDIM1, KDIM2))
      XKER_N_SWETH => PARAM_LIMA_MIXED%XKER_N_SWETH
  CASE('XKER_N_GWETH')
      ALLOCATE(PARAM_LIMA_MIXED%XKER_N_GWETH(KDIM1, KDIM2))
      XKER_N_GWETH => PARAM_LIMA_MIXED%XKER_N_GWETH
  END SELECT
END SUBROUTINE PARAM_LIMA_MIXED_ALLOCATE
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_MIXED
