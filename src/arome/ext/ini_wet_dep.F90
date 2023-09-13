!     ######spl
      SUBROUTINE INI_WET_DEP
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ###########################################################
!
!!****  *INI_RAIN_ICE * - initialize the constants necessary for the warm and
!!                        cold microphysical schemes.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used to
!!    resolve the mixed phase microphysical scheme. The collection kernels of
!!    the precipitating particles are recomputed if necessary if some parameters
!!    defining the ice categories have been modified. The number of small
!!    time steps leading to stable scheme for the rain, ice, snow and ggraupeln
!!    sedimentation is also computed (time-splitting technique).
!!
!!**  METHOD
!!    ------
!!      The constants are initialized to their numerical values and the number
!!    of small time step is computed by dividing the 2* Deltat time interval of
!!    the Leap-frog scheme so that the stability criterion for the rain
!!    sedimentation is fulfilled for a Raindrop maximal fall velocity equal
!!    VTRMAX. The parameters defining the collection kernels are read and are
!!    checked against the new ones. If any change occurs, these kernels are
!!    recomputed and their numerical values are written in the output listiing.
!!
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XP00                 ! Reference pressure
!!        XRD                  ! Gaz constant for dry air
!!        XRHOLW               ! Liquid water density
!!      Module MODD_REF
!!        XTHVREFZ             ! Reference virtual pot.temp. without orography
!!      Module MODD_PARAMETERS_DEP
!!        JPVEXT               !
!!      Module MODD_WET_DEP_DESCR
!!      Module MODD_WET_DEP_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_RAIN_ICE )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!      01-02-2011 M. Mokhtari Adaptation of ini_rain_ice under ini_wet_dep for Aladin
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS_DEP
USE MODD_WET_DEP_DESCR
USE MODD_WET_DEP_PARAM
!
USE MODI_GAMMA
USE MODI_GAMMA_INC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!INTEGER,                 INTENT(IN) :: KLUOUT   ! Logical unit number for prints
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB                ! Coordinates of the first physical
                              ! points along z
INTEGER :: J1,J2              ! Internal loop indexes
REAL :: ZT                    ! Work variable
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZRATE                 ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL :: ZBOUND                ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
REAL :: ZEGS, ZEGR, ZEHS, ZEHG! Bulk collection efficiencies
!
INTEGER :: IND                ! Number of interval to integrate the kernels
REAL :: ZALPHA, ZNU, ZP       ! Parameters to compute the value of the p_moment
                              ! of the generalized Gamma function
REAL :: ZESR                  ! Mean efficiency of rain-aggregate collection
REAL :: ZFDINFTY              ! Factor used to define the "infinite" diameter
!
INTEGER  :: IRESP   ! Return code of FM-routines
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
REAL     :: ZCONC_MAX ! Maximal concentration for snow
REAL     :: ZGAMC,ZGAMC2 ! parameters
                    ! involving various moments of the generalized gamma law
REAL     :: ZFACT_NUCL! Amplification factor for the minimal ice concentration
REAL     :: ZXR     ! Value of x_r in N_r = C_r lambda_r ** x_r
!
INTEGER  :: KND
INTEGER  :: KACCLBDAS,KACCLBDAR,KDRYLBDAG,KDRYLBDAS,KDRYLBDAR
INTEGER  :: KWETLBDAS,KWETLBDAG,KWETLBDAH
REAL     :: PALPHAR,PALPHAS,PALPHAG,PALPHAH
REAL     :: PNUR,PNUS,PNUG,PNUH
REAL     :: PBR,PBS,PBG,PBH
REAL     :: PCR,PCS,PCG,PCH
REAL     :: PDR,PDS,PDG,PDH
REAL     :: PESR,PEGS,PEGR,PEHS,PEHG
REAL     :: PFDINFTY
REAL     :: PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN
REAL     :: PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN
REAL     :: PDRYLBDAR_MAX,PDRYLBDAR_MIN
REAL     :: PWETLBDAS_MAX,PWETLBDAG_MAX,PWETLBDAS_MIN,PWETLBDAG_MIN
REAL     :: PWETLBDAH_MAX,PWETLBDAH_MIN
REAL     :: ZTHVREFZ
!-------------------------------------------------------------------------------
!
!        1.     INTIALIZE OUTPUT LISTING AND COMPUTE KSPLTR FOR EACH MODEL
!               ---------------------------------------------------------
!
!
!*       1.1    Set the hailstones maximum fall velocity
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_WET_DEP',0,ZHOOK_HANDLE)
IF (ALLOCATED(XRTMIN)) THEN       ! In case of nesting microphysics constants of
                                  ! MODD_RAIN_ICE_PARAM are computed only once,
                                  ! but if INI_RAIN_ICE has been called already
                                  ! one must change the XRTMIN size.
  DEALLOCATE(XRTMIN)
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     CHARACTERISTICS OF THE SPECIES
!               ------------------------------
!
!
!*       2.1    Cloud droplet and Raindrop characteristics
!
XAC = (XPI/6.0)*XRHOLW
XBC = 3.0
XCC = XRHOLW*XG/(18.0*1.7E-5) ! Stokes flow (Pruppacher p 322 for T=273K)
XDC = 2.0
!
!
XAR = (XPI/6.0)*XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
!
!XCCR = 1.E7    ! N0_r =  XCXR * lambda_r ** ZXR
XCCR = 8.E6    ! N0_r =  XCXR * lambda_r ** ZXR
ZXR  = -1.     !
!
XF0R = 1.00
XF1R = 0.26
!
XC1R = 1./2.
!
!
!*       2.2    Ice crystal characteristics
!*       2.3    Snowflakes/aggregates characteristics
!
!
XAS = 0.02
XBS = 1.9
XCS = 5.1
XDS = 0.27
!
XCCS = 5.0
XCXS = 1.0
!
XF0S = 0.86
XF1S = 0.28
!
XC1S = 1./XPI
!
!
!*       2.4    Graupel/Frozen drop characteristics
!
!
XAG = 19.6  ! Lump graupel case
XBG = 2.8   ! Lump graupel case
XCG = 124.  ! Lump graupel case
XDG = 0.66  ! Lump graupel case
!
XCCG = 5.E5
XCXG = -0.5
! XCCG = 4.E4 ! Test of Ziegler (1988)
! XCXG = -1.0 ! Test of Ziegler (1988)
!
XF0G = 0.86
XF1G = 0.28
!
XC1G = 1./2.
!
!
!*       2.5    Hailstone characteristics
!
!
XAH = 470.
XBH = 3.0
XCH = 207.
XDH = 0.64
!
!XCCH = 5.E-4
!XCXH = 2.0
!!!!!!!!!!!!
   XCCH = 4.E4 ! Test of Ziegler (1988)
   XCXH = -1.0 ! Test of Ziegler (1988)
!!!    XCCH = 5.E5 ! Graupel_like
!!!    XCXH = -0.5 ! Graupel_like
!!!!!!!!!!!!
!
XF0H = 0.86
XF1H = 0.28
!
XC1H = 1./2.
!
!-------------------------------------------------------------------------------
!
!*       3.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!               ----------------------------------------
!
!
!        3.1    Cloud droplet distribution
!
! Over land
XALPHAC = 1.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
XNUC    = 3.0  ! Gamma law with little dispersion
!
!
! Over sea
XALPHAC2 = 3.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
XNUC2    = 1.0  ! Gamma law with little dispersion
!
!*       3.2    Raindrops distribution
!
XALPHAR = 1.0  ! Exponential law
XNUR    = 1.0  ! Exponential law
!
!*       3.3    Ice crystal distribution
!
XALPHAI = 3.0  ! Gamma law for the ice crystal volume
XNUI    = 3.0  ! Gamma law with little dispersion
!
XALPHAS = 1.0  ! Exponential law
XNUS    = 1.0  ! Exponential law
!
XALPHAG = 1.0  ! Exponential law
XNUG    = 1.0  ! Exponential law
!
XALPHAH = 1.0  ! Gamma law
XNUH    = 8.0  ! Gamma law with little dispersion
!
!*       3.4    Constants for shape parameter
!
ZGAMC = MOMG(XALPHAC,XNUC,3.)
ZGAMC2 = MOMG(XALPHAC2,XNUC2,3.)
XLBC(1)   = XAR*ZGAMC
XLBC(2)   = XAR*ZGAMC2
XLBEXC = 1.0/XBC
!
XLBEXR = 1.0/(-1.0-XBR)
XLBR   = ( XAR*XCCR*MOMG(XALPHAR,XNUR,XBR) )**(-XLBEXR)
!
!XLBEXI = 1.0/(-XBI)
!XLBI   = ( XAI*MOMG(XALPHAI,XNUI,XBI) )**(-XLBEXI)
!
XLBEXS = 1.0/(XCXS-XBS)
XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
!
XLBEXG = 1.0/(XCXG-XBG)
XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG) )**(-XLBEXG)
!
XLBEXH = 1.0/(XCXH-XBH)
XLBH   = ( XAH*XCCH*MOMG(XALPHAH,XNUH,XBH) )**(-XLBEXH)
!
!*       3.5    Minimal values allowed for the mixing ratios
!
XLBDAR_MAX = 100000.0
XLBDAS_MAX = 100000.0
XLBDAG_MAX = 100000.0
!
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
XLBDAS_MAX = ( ZCONC_MAX/XCCS )**(1./XCXS)
!
ALLOCATE( XRTMIN(6) )
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-20
XRTMIN(4) = 1.0E-20
XRTMIN(5) = 1.0E-15
XRTMIN(6) = 1.0E-15
!
XCONC_SEA=1E8 ! 100/cm3
XCONC_LAND=3E8 ! 300/cm3
XCONC_URBAN=5E8 ! 500/cm3
!
!-------------------------------------------------------------------------------
!
!*       4.     CONSTANTS FOR THE SEDIMENTATION
!               -------------------------------
!
!
!*       4.1    Exponent of the fall-speed air density correction
!
XCEXVT = 0.4
!
ZTHVREFZ=300.
ZRHO00 = XP00/(XRD*300.)
!
!*       4.2    Constants for sedimentation
!
XFSEDC(1)  = GAMMA(XNUC+(XDC+3.)/XALPHAC)/GAMMA(XNUC+3./XALPHAC)*     &
            (ZRHO00)**XCEXVT
XFSEDC(2)  = GAMMA(XNUC2+(XDC+3.)/XALPHAC2)/GAMMA(XNUC2+3./XALPHAC2)*     &
            (ZRHO00)**XCEXVT
!
XEXSEDR = (XBR+XDR+1.0)/(XBR+1.0)
XFSEDR  = XCR*XAR*XCCR*MOMG(XALPHAR,XNUR,XBR+XDR)*                         &
            (XAR*XCCR*MOMG(XALPHAR,XNUR,XBR))**(-XEXSEDR)*(ZRHO00)**XCEXVT
!
!  Computations made for Columns
!
XEXRSEDI = 1.9324
XEXCSEDI =-0.9324
XFSEDI   = 3.89745E11*MOMG(XALPHAI,XNUI,3.285)*                          &
                      MOMG(XALPHAI,XNUI,1.7)**(-XEXRSEDI)*(ZRHO00)**XCEXVT
XEXCSEDI =-0.9324*3.0
!
!
XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)
XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
            (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT
!
XEXSEDG = (XBG+XDG-XCXG)/(XBG-XCXG)
XFSEDG  = XCG*XAG*XCCG*MOMG(XALPHAG,XNUG,XBG+XDG)*                         &
            (XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XEXSEDG)*(ZRHO00)**XCEXVT
!
XEXSEDH = (XBH+XDH-XCXH)/(XBH-XCXH)
XFSEDH  = XCH*XAH*XCCH*MOMG(XALPHAH,XNUH,XBH+XDH)*                         &
            (XAH*XCCH*MOMG(XALPHAH,XNUH,XBH))**(-XEXSEDH)*(ZRHO00)**XCEXVT
!
!
!-------------------------------------------------------------------------------
!
!*       5.     CONSTANTS FOR THE SLOW COLD PROCESSES
!               -------------------------------------
!
!
!*       5.1    Constants for ice nucleation
!
!
!
XMNU0 = 6.88E-13
!
!*       5.2    Constants for vapor deposition on ice
!
XSCFAC = (0.63**(1./3.))*SQRT((ZRHO00)**XCEXVT) ! One assumes Sc=0.63
!
X0DEPI = (4.0*XPI)*XC1I*XF0I*MOMG(XALPHAI,XNUI,1.)
X2DEPI = (4.0*XPI)*XC1I*XF2I*XC_I*MOMG(XALPHAI,XNUI,XDI+2.0)
!
X0DEPS = (4.0*XPI)*XCCS*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0*XPI)*XCCS*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = XCXS-1.0
XEX1DEPS = XCXS-0.5*(XDS+3.0)
!
X0DEPG = (4.0*XPI)*XCCG*XC1G*XF0G*MOMG(XALPHAG,XNUG,1.)
X1DEPG = (4.0*XPI)*XCCG*XC1G*XF1G*SQRT(XCG)*MOMG(XALPHAG,XNUG,0.5*XDG+1.5)
XEX0DEPG = XCXG-1.0
XEX1DEPG = XCXG-0.5*(XDG+3.0)
!
X0DEPH = (4.0*XPI)*XCCH*XC1H*XF0H*MOMG(XALPHAH,XNUH,1.)
X1DEPH = (4.0*XPI)*XCCH*XC1H*XF1H*SQRT(XCH)*MOMG(XALPHAH,XNUH,0.5*XDH+1.5)
XEX0DEPH = XCXH-1.0
XEX1DEPH = XCXH-0.5*(XDH+3.0)
!
!*       5.3    Constants for pristine ice autoconversion
!
XTIMAUTI = 1.E-3  !  Time constant at T=T_t
XTEXAUTI = 0.015  !  Temperature factor of the I+I collection efficiency
!!XCRIAUTI = 0.25E-3 !  Critical ice content for the autoconversion to occur
XCRIAUTI = 0.2E-4 !  Critical ice content for the autoconversion to occur
                  !  Revised value by Chaboureau et al. (2001)
!
!
!*       5.4    Constants for snow aggregation
!
XCOLIS   = 0.25 ! Collection efficiency of I+S
XCOLEXIS = 0.05 ! Temperature factor of the I+S collection efficiency
XFIAGGS  = (XPI/4.0)*XCOLIS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = XCXS-XDS-2.0
!
!-------------------------------------------------------------------------------
!
!*       6.     CONSTANTS FOR THE SLOW WARM PROCESSES
!               -------------------------------------
!
!
!*       6.1    Constants for the cloud droplets autoconversion
!
XTIMAUTC = 1.E-3
XCRIAUTC = 0.5E-3
!
!*       6.2    Constants for the accretion of cloud droplets by raindrops
!
XFCACCR  = (XPI/4.0)*XCCR*XCR*(ZRHO00**XCEXVT)*MOMG(XALPHAR,XNUR,XDR+2.0)
XEXCACCR = -XDR-3.0
!
!*       6.3    Constants for the evaporation of the raindrops
!
X0EVAR = (4.0*XPI)*XCCR*XC1R*XF0R*MOMG(XALPHAR,XNUR,1.)
X1EVAR = (4.0*XPI)*XCCR*XC1R*XF1R*SQRT(XCR)*MOMG(XALPHAR,XNUR,0.5*XDR+1.5)
XEX0EVAR = -2.0
XEX1EVAR = -1.0-0.5*(XDR+3.0)
!
!
IF (LHOOK) CALL DR_HOOK('INI_WET_DEP',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
CONTAINS

!
!*       11.    INTERNAL FUNCTIONS
!               -------------------
!
!
!*       11.1   p_moment of the Generalized GAMMA function
!
REAL FUNCTION MOMG(ZALPHA,ZNU,ZP)

  IMPLICIT NONE
  REAL, INTENT(IN) :: ZALPHA,ZNU,ZP

  MOMG = GAMMA(ZNU+ZP/ZALPHA)/GAMMA(ZNU)

END FUNCTION MOMG
!
!


END SUBROUTINE INI_WET_DEP
