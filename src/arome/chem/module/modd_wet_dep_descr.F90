!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ##########################
      MODULE MODD_WET_DEP_DESCR
!     ##########################
!
!!****  *MODD_WET_DEP_DESCR* - declaration of the microphysical descriptive
!!                              constants for use in the warm and cold schemes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop and
!     the ice crystal habits and the parameters relevant of the dimensional
!     distributions.
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
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_RAIN_ICE_DESCR)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!       J.-P. Pinty   29/11/02 add ICE4
!!      01-02-2011 M. Mokhtari Adaptation of MODD_RAIN_ICE_DESCR under MODD_WET_DEP_DESCR for Aladin
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
REAL,SAVE :: XCEXVT               ! air density fall speed correction
!
REAL,SAVE :: XAC,XBC,XCC,XDC                          ! Cloud droplet  charact.
REAL,SAVE :: XAR,XBR,XCR,XDR,XCCR     ,XF0R,XF1R,XC1R ! Raindrop       charact.
REAL,SAVE :: XAI,XBI,XC_I,XDI          ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL,SAVE :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
REAL,SAVE :: XAG,XBG,XCG,XDG,XCCG,XCXG,XF0G,XF1G,XC1G ! Graupel        charact.
REAL,SAVE :: XAH,XBH,XCH,XDH,XCCH,XCXH,XF0H,XF1H,XC1H ! Hail           charact.
!
REAL,SAVE :: XALPHAC,XNUC,XALPHAC2,XNUC2, XLBEXC      ! Cloud droplet  distribution parameters
REAL,DIMENSION(2), SAVE :: XLBC ! Cloud droplet distribution parameters
REAL,SAVE :: XALPHAR,XNUR,XLBEXR,XLBR ! Raindrop       distribution parameters
REAL,SAVE :: XALPHAI,XNUI,XLBEXI,XLBI ! Cloud ice      distribution parameters
REAL,SAVE :: XALPHAS,XNUS,XLBEXS,XLBS ! Snow/agg.      distribution parameters
REAL,SAVE :: XALPHAG,XNUG,XLBEXG,XLBG ! Graupel        distribution parameters
REAL,SAVE :: XALPHAH,XNUH,XLBEXH,XLBH ! Hail           distribution parameters
!
REAL,SAVE :: XLBDAR_MAX,XLBDAS_MAX,XLBDAG_MAX ! Max values allowed for the shape
                                              ! parameters (rain,snow,graupeln)
!
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XRTMIN ! Min values allowed for the mixing ratios
REAL,DIMENSION(:,:,:), SAVE, ALLOCATABLE :: XCONC ! Concentration of cloud droplet
REAL,SAVE :: XCONC_SEA   ! Diagnostic concentration of droplets over sea
REAL,SAVE :: XCONC_LAND  !  Diagnostic concentration of droplets over land
REAL,SAVE :: XCONC_URBAN ! Diagnostic concentration of droplets over urban area
!
END MODULE MODD_WET_DEP_DESCR
