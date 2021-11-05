!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_INI_CST
!     ###################
!
INTERFACE
!
SUBROUTINE INI_CST
END SUBROUTINE INI_CST 
!
END INTERFACE
!
END MODULE MODI_INI_CST
!
!
!
!     ##################
      SUBROUTINE INI_CST 
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
!!      J. Escobar  28/03/2014 for pb with emissivity/aerosol reset XMNH_TINY=1.0e-80 in real8 case 
!!      J.Escobar : 10/2017 : for real*4 , add XMNH_HUGE_12_LOG
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      J.Escobar : 5/10/2018 : for real*4 ,higher value for XEPS_DT = 1.5e-4
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
use modd_precision, only: MNHREAL
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*	 1.     FUNDAMENTAL CONSTANTS
!	        ---------------------
!
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
!JUAN OVERFLOW XSTEFAN = 2.* XPI**5 * XBOLTZ**4 / (15.* XLIGHTSPEED**2 * XPLANCK**3)
XSTEFAN = ( 2.* XPI**5 / 15. ) * ( (XBOLTZ / XPLANCK) * XBOLTZ ) * (XBOLTZ/(XLIGHTSPEED*XPLANCK))**2 
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
!   Some machine precision value depending of real4/8 use  
!


XMNH_EPSILON = EPSILON (XMNH_EPSILON )
XMNH_HUGE    = HUGE    (XMNH_HUGE )
XMNH_HUGE_12_LOG = LOG ( SQRT(XMNH_HUGE)  )

#if (MNH_REAL == 8)
XMNH_TINY      = 1.0e-80_MNHREAL
XEPS_DT        = 1.0e-5_MNHREAL
XRES_FLAT_CART = 1.0e-12_MNHREAL
XRES_OTHER     = 1.0e-9_MNHREAL
XRES_PREP      = 1.0e-8_MNHREAL
#elif (MNH_REAL == 4)
XMNH_TINY      = TINY    (XMNH_TINY    )
XEPS_DT        = 1.5e-4_MNHREAL
XRES_FLAT_CART = 1.0e-12_MNHREAL
XRES_OTHER     = 1.0e-7_MNHREAL
XRES_PREP      = 1.0e-4_MNHREAL
#else
#error "Invalid MNH_REAL"
#endif
XMNH_TINY_12 = SQRT    (XMNH_TINY    )



!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_CST 
