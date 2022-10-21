!     ############################
      MODULE MODD_PARAM_LIMA_MIXED
!     ###########################{
!
!!****  *MODD_PARAM_LIMA_MIXED* - declaration of some descriptive parameters and
!!                                microphysical factors extensively used in 
!!                                the LIMA mixed scheme.
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
REAL,SAVE :: XAG,XBG,XCG,XDG,XCCG,XCXG,XF0G,XF1G,XC1G ! Graupel        charact.
REAL,SAVE :: XLBEXG,XLBG              ! Graupel        distribution parameters 
REAL,SAVE :: XLBDAG_MAX               ! Max values allowed for the shape
                                      ! parameter of graupeln
!
REAL,SAVE :: XAH,XBH,XCH,XDH,XCCH,XCXH,XF0H,XF1H,XC1H ! Hail           charact.
REAL,SAVE :: XALPHAH,XNUH,XLBEXH,XLBH ! Hail           distribution parameters
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS - Graupel
!             -------------------------------
!
REAL,SAVE :: XFSEDG, XEXSEDG                     ! Sedimentation fluxes of Graupel
!
REAL,SAVE :: X0DEPG,X1DEPG,XEX0DEPG,XEX1DEPG     ! Deposition on graupel
!
REAL,SAVE :: XHMTMIN,XHMTMAX,XHM1,XHM2,        & ! Constants for the
             XHM_YIELD,XHM_COLLCS,XHM_FACTS,   & !      revised
                       XHM_COLLCG,XHM_FACTG,   & ! Hallett-Mossop process
    	     XGAMINC_HMC_BOUND_MIN,            & ! Min val. of Lbda_c for HMC
    	     XGAMINC_HMC_BOUND_MAX,            & ! Max val. of Lbda_c for HMC
             XHMSINTP1,XHMSINTP2,              & ! (this is no more used !)
             XHMLINTP1,XHMLINTP2
!
REAL,SAVE :: XDCSLIM,XCOLCS,                   & ! Constants for the riming of
    	     XEXCRIMSS,XCRIMSS,                & ! the aggregates : RIM
    	     XEXCRIMSG,XCRIMSG,                & !
    	     XEXSRIMCG,XSRIMCG,                & !
    	     XGAMINC_BOUND_MIN,                & ! Min val. of Lbda_s for RIM
    	     XGAMINC_BOUND_MAX,                & ! Max val. of Lbda_s for RIM
    	     XRIMINTP1,XRIMINTP2                 ! Csts for lin. interpol. of 
                                                 ! the tab. incomplete Gamma law
INTEGER,SAVE :: NGAMINC                          ! Number of tab. Lbda_s
REAL, DIMENSION(:), SAVE, ALLOCATABLE          &
                       :: XGAMINC_RIM1,        & ! Tab. incomplete Gamma funct.
                          XGAMINC_RIM2,        & ! for XDS+2 and for XBS
                          XGAMINC_HMC            ! and for the HM process
!
REAL,SAVE :: XFRACCSS,                         & ! Constants for the accretion 
             XLBRACCS1,XLBRACCS2,XLBRACCS3,    & ! raindrops onto the aggregates
             XFSACCRG,                         & ! ACC (processes RACCSS and
             XLBSACCR1,XLBSACCR2,XLBSACCR3,    & !                SACCRG)
    	     XACCLBDAS_MIN,                    & ! Min val. of Lbda_s for ACC
    	     XACCLBDAS_MAX,                    & ! Max val. of Lbda_s for ACC
    	     XACCLBDAR_MIN,                    & ! Min val. of Lbda_r for ACC
    	     XACCLBDAR_MAX,                    & ! Max val. of Lbda_r for ACC
    	     XACCINTP1S,XACCINTP2S,            & ! Csts for bilin. interpol. of 
    	     XACCINTP1R,XACCINTP2R               !   Lbda_s and Lbda_r in the 
                        			 ! XKER_RACCSS and XKER_SACCRG
                          			 !            tables
INTEGER,SAVE :: NACCLBDAS,                     & ! Number of Lbda_s values and
    	        NACCLBDAR                        !   of Lbda_r values in the
                        			 ! XKER_RACCSS and XKER_SACCRG
                        			 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
       			 :: XKER_RACCSS,       & ! Normalized kernel for RACCSS
       			    XKER_RACCS,        & ! Normalized kernel for RACCS
       			    XKER_SACCRG          ! Normalized kernel for SACCRG 
REAL,SAVE :: XFSCVMG                             ! Melting-conversion factor of
                                                 ! the aggregates
!
REAL,SAVE :: XCOLIR,                           & ! Constants for rain contact
    	     XEXRCFRI,XRCFRI,                  & ! freezing : CFR
    	     XEXICFRR,XICFRR                     !
!
REAL,SAVE :: XFCDRYG,                          & ! Constants for the dry growth
             XCOLCG,                           & ! of the graupeln :
             XCOLIG,XCOLEXIG,XFIDRYG,          & ! 
             XCOLSG,XCOLEXSG,XFSDRYG,          & !             RCDRYG
             XLBSDRYG1,XLBSDRYG2,XLBSDRYG3,    & !             RIDRYG
             XFRDRYG,                          & !             RSDRYG
             XLBRDRYG1,XLBRDRYG2,XLBRDRYG3,    & !             RRDRYG
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
INTEGER,SAVE :: NDRYLBDAR,                     & ! Number of Lbda_r,
    	        NDRYLBDAS,                     & !        of Lbda_s and
    	        NDRYLBDAG                        !        of Lbda_g values in
    	                                         ! the XKER_SDRYG and XKER_RDRYG
                        			 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_SDRYG,        & ! Normalized kernel for SDRYG
       			    XKER_RDRYG           ! Normalized kernel for RDRYG 
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS - Hail
!             ----------------------------
!
REAL,SAVE :: XFSEDH,XEXSEDH                      ! Constants for sedimentation
!
!
REAL,SAVE :: X0DEPH,X1DEPH,XEX0DEPH,XEX1DEPH     ! Constants for deposition
!
REAL,SAVE :: XFWETH,XFSWETH,                   & ! Constants for the wet growth
             XLBSWETH1,XLBSWETH2,XLBSWETH3,    & ! of the hailstones : WET
             XFGWETH,                          & !   processes RSWETH
             XLBGWETH1,XLBGWETH2,XLBGWETH3,    & !             RGWETH
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
INTEGER,SAVE :: NWETLBDAS,                     & ! Number of Lbda_s,
                NWETLBDAG,                     & !        of Lbda_g and
                NWETLBDAH                        !        of Lbda_h values in
                                                 ! the XKER_SWETH and XKER_GWETH
                                                 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_SWETH,        & ! Normalized kernel for SWETH
                            XKER_GWETH           ! Normalized kernel for GWETH

!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_MIXED
