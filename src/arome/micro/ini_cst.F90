!     ######spl
      SUBROUTINE INI_CST 
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################
!
!!****  *INI_CST * - routine to initialize the module MODD_CST
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the physical constants
!     stored in  module MODD_CST.
!      
!
!!**  METHOD
!!    ------
!!      The physical constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!      FMLOOK : to retrieve logical unit number associated to a file
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST     : contains physical constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CST, routine INI_CST)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/05/94 
!!      J. Stein    02/01/95  add the volumic mass of liquid water
!!      J.-P. Pinty 13/12/95  add the water vapor pressure over solid ice
!!      J. Stein    29/06/97  add XTH00
!!      V. Masson   05/10/98  add XRHOLI
!!      C. Mari     31/10/00  add NDAYSEC
!!      V. Masson   01/03/03  add XCONDI
!!      R. El Khatib 04/08/14 add pre-computed quantities
!!      P. Marguinaud 04/10/16 Port to single precision
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*	 1.     FUNDAMENTAL CONSTANTS
!	        ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_CST',0,ZHOOK_HANDLE)
XPI         = 2.*ASIN(1.)
XKARMAN     = 0.4
XLIGHTSPEED = 299792458.
XPLANCK     = 6.6260755E-34
XBOLTZ      = 1.380658E-23
XAVOGADRO   = 6.0221367E+23
!
!-------------------------------------------------------------------------------
!
!*       2.     ASTRONOMICAL CONSTANTS
!	        ----------------------
!
XDAY   = 86400.
XSIYEA = 365.25*XDAY*2.*XPI/ 6.283076
XSIDAY = XDAY/(1.+XDAY/XSIYEA)
XOMEGA = 2.*XPI/XSIDAY
NDAYSEC = 24*3600 ! Number of seconds in a day
!
!-------------------------------------------------------------------------------!
!
!
!*       3.     TERRESTRIAL GEOIDE CONSTANTS
!	        ----------------------------
!
XRADIUS = 6371229.
XG      = 9.80665
!
!-------------------------------------------------------------------------------
!
!*	 4.     REFERENCE PRESSURE
!	        -------------------
!
! Ocean model cst same as in 1D/CMO SURFEX
! values used in ini_cst to overwrite XP00 and XTH00
XRH00OCEAN =1024.
XTH00OCEAN = 286.65
XSA00OCEAN= 32.6
XP00OCEAN = 201.E5
!Atmospheric model
XP00 = 1.E5
XTH00 = 300.
!-------------------------------------------------------------------------------
!
!*	 5.     RADIATION CONSTANTS
!	        -------------------
!
XSTEFAN = REAL (2._8* REAL (XPI, 8)**5 * REAL (XBOLTZ, 8)**4 / &
        & (15._8* REAL (XLIGHTSPEED, 8)**2 * REAL (XPLANCK, 8)**3))
XI0     = 1370.
!
!-------------------------------------------------------------------------------
!
!*	 6.     THERMODYNAMIC CONSTANTS
!	        -----------------------
!
XMD    = 28.9644E-3
XMV    = 18.0153E-3
XRD    = XAVOGADRO * XBOLTZ / XMD
XRV    = XAVOGADRO * XBOLTZ / XMV
XEPSILO= XMV/XMD
XCPD   = 7.* XRD /2.
XCPV   = 4.* XRV
XRHOLW = 1000.
XRHOLI = 900.
XCONDI = 2.22
XCL    = 4.218E+3
XCI    = 2.106E+3
XTT    = 273.16
XLVTT  = 2.5008E+6
XLSTT  = 2.8345E+6
XLMTT  = XLSTT - XLVTT
XESTT  = 611.14
XGAMW  = (XCL - XCPV) / XRV
XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
XALPW  = LOG(XESTT) + (XBETAW /XTT) + (XGAMW *LOG(XTT))
XGAMI  = (XCI - XCPV) / XRV
XBETAI = (XLSTT/XRV) + (XGAMI * XTT)
XALPI  = LOG(XESTT) + (XBETAI /XTT) + (XGAMI *LOG(XTT))
! Values identical to ones used in CMO1D in SURFEX /could be modified
! Coefficient of thermal expansion of water (K-1)
XALPHAOC = 1.9E-4
! Coeff of Haline contraction coeff (S-1)
XBETAOC= 7.7475E-4
!
!*	 7.     PRECOMPUTED CONSTANTS
!	        ---------------------
!
RDSRV = XRD/XRV
RDSCPD = XRD/XCPD
RINVXP00 =  1./XP00
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INI_CST',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CST 
